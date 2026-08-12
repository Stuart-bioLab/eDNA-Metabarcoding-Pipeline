import argparse
import configparser
from pathlib import Path
from datetime import datetime
import logging
import subprocess
from Bio.Seq import Seq
import sys
import shutil

def load_config(config_file):
    """Load configuration file."""
    config = configparser.ConfigParser()
    config.read(config_file) # read in default arguments from config file
    return config

def get_args(config):
    """Read in command line arguments and set config file defaults as values if none are given."""
    parser = argparse.ArgumentParser(
        description="Run eDNA pipeline"
    )

    parser.add_argument(
        "--manifest",
        default=config["PATHS"]["manifest"],
        help="manifest file mapping sample ids to read paths"
    )

    parser.add_argument(
        "--metadata",
        default=config["PATHS"]["metadata"],
        help="metadata for samples"
    )

    parser.add_argument(
        "--threads",
        default=config["PARAMETERS"]["threads"],
        help="Number of CPU cores to use for applicable functions"
    )

    parser.add_argument( # MIGHT ONLY KEEP THIS FOR DEVELOPMENT
        "--archive",
        help="supply archive file to skip importing step"
    )

    parser.add_argument(
        "--forward_primer",
        default=config["PRIMERS"]["forward"],
        help="forward primer sequence"
    )

    parser.add_argument(
        "--reverse_primer",
        default=config["PRIMERS"]["reverse"],
        help="reverse primer sequence"
    )

    parser.add_argument(
        "--trim_forward",
        default=config["PARAMETERS"]["trim_forward"],
        help="quality info for dada2. where to trim forward reads before denoising"
    )

    parser.add_argument(
        "--trim_reverse",
        default=config["PARAMETERS"]["trim_reverse"],
        help="quality info for dada2. where to trim reverse reads before denoising"
    )

    parser.add_argument(
        "--trunc_forward",
        default=config["PARAMETERS"]["trunc_forward"],
        help="quality info for dada2. where to truncate forward reads before denoising"
    )

    parser.add_argument(
        "--trunc_reverse",
        default=config["PARAMETERS"]["trunc_reverse"],
        help="quality info for dada2. where to truncate reverse reads before denoising"
    )
    
    parser.add_argument(
        "--database_tax",
        default=config["PATHS"]["database_tax"],
        help="refrence database taxa file in QIIME format for vsearch and naive bayes taxonomy assignment steps"
    )
    
    parser.add_argument(
        "--database_seq",
        default=config["PATHS"]["database_seq"],
        help="refrence database seq file in QIIME format for vsearch and naive bayes taxonomy assignment steps"
    )

    parser.add_argument( # YOU CAN GET RID OF THIS LATER
        "--asv_seqs",
        help="for development: skip up to taxonomy assignment if denoised seqs provided"
    )

    # parser.add_argument("-bdt", "--blast_database_tax", help="BLAST database taxa file in QIIME format for BLAST taxonomy assignment")
    # parser.add_argument("-bds", "--blast_database_seq", help="BLAST database deq file in QIIME format for BLAST taxonomy assignment")

    return parser.parse_args()

def create_outdir(base):
    """Create output directory named with the current time to differentiate runs."""
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(base) / run_id
    outdir.mkdir(parents=True) # make results dir if it does not exist yet
    return outdir

def setup_logger(log_file):
    """Initialize log file."""
    logger = logging.getLogger(__name__) # initialize logger and name it after this script
    logger.setLevel(logging.INFO) # show everything but debug information

    file_handler = logging.FileHandler(log_file) # send log messages to this file
    file_handler.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter) # show time, level, and message in the log

    logger.addHandler(file_handler) # tell logger where to send messages

    return logger

def import_reads(logger, manifest, archive):
    """Import reads into qiime2 from manifest."""
    logger.info("importing reads")

    try:
        subprocess.run(
            [
                "qiime", "tools", "import",
                "--type", "SampleData[PairedEndSequencesWithQuality]",
                "--input-path", manifest,
                "--input-format", "PairedEndFastqManifestPhred33V2",
                "--output-path", archive
            ],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error("import failed")
        logger.error(e.stderr)
        sys.exit(1)
    
    logger.info("DONE importing reads")
    return

def trim_reads(logger, primers, reads_archive, threads, outdir):
    forward_primer, reverse_primer = primers

    rev_comp_forward = str(Seq(forward_primer).reverse_complement())
    rev_comp_reverse = str(Seq(reverse_primer).reverse_complement())

    logger.info(f"loaded primers: {forward_primer}, {reverse_primer}")
    logger.info(f"loaded reverse comp primers: {rev_comp_forward}, {rev_comp_reverse}")

    reverse_trimmed_reads = outdir / "reverse_trimmed_reads.qza"
    trimmed_reads = outdir / "trimmed_reads.qza"

    logger.info("reverse trimming reads")

    try:
        subprocess.run(
            [
                "qiime", "cutadapt", "trim-paired",
                "--i-demultiplexed-sequences", reads_archive,
                "--p-adapter-f", rev_comp_forward,
                "--p-adapter-r", rev_comp_reverse,
                "--p-match-read-wildcards", "TRUE",
                "--p-error-rate", "0.25", 
                "--p-cores", threads,
                "--o-trimmed-sequences", reverse_trimmed_reads,
                "--o-stats", outdir / "reverse_trim_stats.qza"
            ],
            capture_output=True,
            text=True,
            check=True
        )
    
    except subprocess.CalledProcessError as e:
        logger.error("reverse trimming failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info("DONE reverse trimming reads")

    logger.info("trimming reads")

    try:
        subprocess.run(
            [
                "qiime", "cutadapt", "trim-paired",
                "--i-demultiplexed-sequences", reverse_trimmed_reads,
                "--p-front-f", forward_primer,
                "--p-front-r", reverse_primer,
                "--p-cores", threads,
                "--p-match-read-wildcards", "TRUE",
                "--p-discard-untrimmed", "TRUE",
                "--p-match-adapter-wildcards", "TRUE",
                "--o-trimmed-sequences", trimmed_reads,
                "--o-stats", outdir / "trim_stats.qza"
            ],
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        logger.error("trimming failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info("DONE trimming reads")

    return trimmed_reads

def denoise_reads(logger, reads, params, outdir):
    """Run DADA2 to denoise reads and generate feature table."""
    asv_seqs = outdir / "asv_seqs.qza"
    feat_table = outdir / "feat_table"
    denoise_stats = outdir / "denoise_stats"
    transition_stats = outdir / "transition_stats"

    logger.info("denoising reads")

    try: 
        subprocess.run(
            [
                "qiime", "dada2", "denoise-paired",
                "--i-demultiplexed-seqs", reads,
                "--p-trim-left-f", params["trim_forward"],
                "--p-trim-left-r", params["trim_reverse"],
                "--p-trunc-len-f", params["trunc_forward"],
                "--p-trunc-len-r", params["trunc_reverse"],
                "--p-n-threads", params["threads"],
                "--o-representative-sequences", asv_seqs,
                "--o-table", feat_table,
                "--o-denoising-stats", denoise_stats,
                "--o-base-transition-stats", transition_stats
            ],
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        logger.error("denoising failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info("DONE denoising reads")

    return asv_seqs, feat_table

def run_vsearch(logger, asv_seqs, ref_tax, ref_seq, threads, outdir):
    """Get exact matches with vsearch"""
    out_tax = outdir / "vsearch_tax.qza"
    top_hits = outdir / "vsearch_top_hits.qza"

    logger.info("running vsearch")

    try:
        subprocess.run(
            [
                "qiime", "feature-classifier", "classify-consensus-vsearch",
                "--i-query", asv_seqs,
                "--i-reference-reads", ref_seq,
                "--i-reference-taxonomy", ref_tax,
                "--p-perc-identity", "1.0",
                "--p-min-consensus", "0.94",
                "--p-threads", threads,
                "--o-classification", out_tax,
                "--o-search-results", top_hits
            ],
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        logger.error("vsearch failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info("DONE running vsearch")

    return out_tax

def extract_qza(logger, archive, new_name, outdir):
    """Unzip qiime archive. Allows for further data manipulation. Much faster than the export command."""
    tmpdir = outdir / "tmp"
    tmpdir.mkdir()

    logger.info(f"extracting {archive}")
    try:
        subprocess.run(
            [
                "unzip",
                "-qd",
                tmpdir,
                archive
            ],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"error: {e}")
        sys.exit(1)

    path_to_data = list(tmpdir.glob("*/data/*"))[0]
    final_path = outdir / new_name
    try:
        subprocess.run(
            [
                "mv",
                path_to_data,
                final_path
            ],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"error: {e}")
        sys.exit(1)
    logger.info(f"DONE extracting {archive}")

    shutil.rmtree(tmpdir)
    return final_path

def main():
    config = load_config("config.ini")
    args = get_args(config)
    outdir = create_outdir("results")
    logger = setup_logger(outdir / "pipeline.log")

    threads = args.threads
    logger.info(f"starting pipeline with {threads} threads")

    metadata = Path(args.metadata).resolve()
    logger.info(f"loaded metadata: {metadata}")

    if not args.archive: # only import reads if archive not supplied (this takes a while)
        reads_archive = Path(outdir / "reads.qza").resolve()
        manifest = Path(args.manifest).resolve() # only load manifest if there's no archive supplied
        logger.info(f"loaded manifest: {manifest}")
        import_reads(logger, manifest, reads_archive)
    else:
        logger.info("skipping import and using archive file instead")
        reads_archive = Path(args.archive).resolve() # get archive path
    logger.info(f"loaded archive: {reads_archive}")

    if not args.asv_seqs:
        primers = args.forward_primer, args.reverse_primer
        trim_dir = outdir / "trimmed_reads" # store trimmed read file here
        trim_dir.mkdir()
        trimmed_reads = trim_reads(logger, primers, reads_archive, threads, trim_dir)

        dada_params = { # dada2 parameters from command line arguments
            "trim_forward": args.trim_forward,
            "trim_reverse": args.trim_reverse,
            "trunc_forward": args.trunc_forward,
            "trunc_reverse": args.trunc_reverse,
            "threads": threads
        }
        dada_dir = outdir / "dada2" # store deniosing files here
        dada_dir.mkdir()
        asv_seqs, feat_table = denoise_reads(logger, trimmed_reads, dada_params, dada_dir)
    else:
        logger.info("skipping denoising and using provided seqs instead")
        asv_seqs = Path(args.asv_seqs).resolve()
        logger.info(f"loaded seqs: {asv_seqs}")

    crabs_ref_tax = args.database_tax
    crabs_ref_seq = args.database_seq
    logger.info("loaded crabs reference db files")
    vsearch_dir = outdir / "vsearch"
    vsearch_dir.mkdir()
    vsearch_out = run_vsearch(logger, asv_seqs, crabs_ref_tax, crabs_ref_seq, threads, vsearch_dir)

    new_name = "vsearch_tax.tsv"
    print(extract_qza(logger, vsearch_out, new_name, vsearch_dir))

    logger.info("pipeline end")

if __name__ == "__main__":
    main()