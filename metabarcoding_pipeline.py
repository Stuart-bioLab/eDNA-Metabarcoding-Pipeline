import argparse
import configparser
from pathlib import Path
from datetime import datetime
import logging
import subprocess
from Bio.Seq import Seq
import sys
import shutil
import pandas as pd

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
        default=None,
        help="manifest file mapping sample ids to read paths"
    )

    parser.add_argument(
        "--threads",
        default=config["PARAMETERS"]["threads"],
        help="Number of CPU cores to use for applicable functions"
    )

    parser.add_argument( # MIGHT ONLY KEEP THIS FOR DEVELOPMENT
        "--archive",
        default=None,
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
        "--crabs_database_tax",
        default=config["DATABASES"]["crabs_database_tax"],
        help="refrence database taxa file in QIIME format for vsearch and naive bayes taxonomy assignment steps"
    )
    
    parser.add_argument(
        "--crabs_database_seq",
        default=config["DATABASES"]["crabs_database_seq"],
        help="refrence database seq file in QIIME format for vsearch and naive bayes taxonomy assignment steps"
    )

    parser.add_argument(
        "--blast_database_tax",
        default=config["DATABASES"]["blast_database_tax"],
        help="refrence database tax file in QIIME format for BLAST"
    )

    parser.add_argument(
        "--blast_database_seq",
        default=config["DATABASES"]["blast_database_seq"],
        help="refrence database seq file in QIIME format for BLAST"
    )

    parser.add_argument(
        "--perc_identity",
        default=config["PARAMETERS"]["perc_identity"],
        help="minimum percent identity accepted for blast hits"
    )

    parser.add_argument(
        "--query_cov",
        default=config["PARAMETERS"]["query_cov"],
        help="minimum alignment coverage accepted for blast hits"
    )

    parser.add_argument(
        "--max_accepts",
        default=config["PARAMETERS"]["max_accepts"],
        help="maximum number of hits to keep for blast query"
    )

    parser.add_argument( # FOR DEV PURPOSES
        "--skip_to_mapping",
        action="store_true"
    )

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
    """Trim primers and reverse complemnt primers."""
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
    feat_table = outdir / "feat_table.qza"
    denoise_stats = outdir / "denoise_stats"

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
    """Get exact matches with VSEASRCH."""
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
            ["unzip", "-qd", tmpdir, archive],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"extracting {archive} failed")
        logger.error(e.stderr)
        sys.exit(1)

    path_to_data = list(tmpdir.glob("*/data/*"))[0]
    final_path = outdir / new_name
    try:
        subprocess.run(
            ["mv", path_to_data, final_path],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"moving {archive} failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info(f"extracted {archive} to {final_path}")

    shutil.rmtree(tmpdir)
    return final_path

def parse_output(logger, tax_file, outdir, level):
    """Parse output from taxa classification and separate into retained and unassigned tsvs."""
    prefix = tax_file.name.split("_")[0] # get the prefix of the input file to name the output

    new_name = prefix + "_tax.tsv"
    unzipped_tax = extract_qza(logger, tax_file, new_name, outdir) # extract taxonomy tsv from qiime archive

    tax_df = pd.read_csv(unzipped_tax, sep="\t") # read in vsearch output taxonomy.tsv

    logger.info(f"parsing {prefix} taxonomy")

    unassigned = tax_df[
        tax_df["Taxon"] # take the column with taxonomic classification
        .str.split(";", expand=False) # split it into a list of taxonomic levels
        .apply(lambda x: len(x) < level) # if any rows have len < level, they are not classified deep enough
    ]
    
    unassigned_out = outdir / (prefix + "_unassigned_tax.tsv") # this goes into the map_seqs function to figure out which exact seqs were left unassigned so they can be sent to the next step
    unassigned.to_csv(unassigned_out, sep="\t", index=False)

    retained = tax_df.drop( # drop all rows shared with unassigned (retain only family-level classifications)
        tax_df[
            tax_df["Feature ID"]
            .isin(unassigned["Feature ID"])
        ]
        .index
    )
    
    retained_out = outdir / (prefix + "_retained_tax.tsv") # this is used to prep taxonomy file for phyloseq
    retained.to_csv(retained_out, sep="\t", index=False)

    logger.info(f"DONE parsing {prefix} taxonomy")

    return unassigned_out, retained_out

def filter_seqs(logger, asv_seqs, unassigned, outdir):
    """Use unassigned ids to find unassigned sequences. Filter for unassigned sequences and re-import for next step."""
    prefix = unassigned.name.split("_")[0]

    new_name = prefix + "_seq.fasta"
    unzipped_seq = extract_qza(logger, asv_seqs, new_name, outdir)

    logger.info(f"filtering {prefix} seqs")

    asv_map = {}
    with open(unzipped_seq) as f: # map asv sequences to their qiime feature id
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            feat = lines[i][1:].strip()
            seq = lines[i+1].strip()
            asv_map[feat] = seq

    features_index = list(pd.read_csv(unassigned, sep="\t")["Feature ID"]) # get all ids that were unassigned
    
    unassigned_fasta = outdir / (prefix + "_unassigned_seq.fasta") # non-archive version of below which i pass to phyloseq processing part to trim taxa that are still unassigned after bayesian classification
    with open(unassigned_fasta, "w") as f: # write out file with only the sequences corresponding to unassigned features
        for feat in features_index:
            f.write(f">{feat}\n")
            f.write(f"{asv_map[feat]}\n")

    logger.info(f"DONE filtering {prefix} sequences")
    logger.info(f"importing filtered {prefix} seqs")

    asv_out =  outdir / (prefix + "_unassigned_seq.qza") # unassigned seqs that go to the next classifier 

    try:
        subprocess.run( # import filtered sequences
            [
                "qiime", "tools", "import",
                "--input-path", unassigned_fasta,
                "--output-path", asv_out,
                "--type", "FeatureData[Sequence]"
            ],
            capture_output=True,
            text=True,
            check=True
        )
    
    except subprocess.CalledProcessError as e:
        logger.error("importing failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info(f"DONE importing filtered {prefix} seqs")

    return asv_out, unassigned_fasta

def run_nb_classifier(logger, asv_seqs, train_tax, train_seq, threads, outdir):
    """Train and run the naive Bayes classifier on sequences vsearch did not have exact matches for."""
    rescript_classifier = outdir / "rescript_classifier"
    rescript_evaluation = outdir / "rescript_evaluation"
    rescript_observed_taxonomy = outdir / "rescript_observed_taxonomy"

    logger.info("training nb classifier")
    try:
        subprocess.run(
            [
                "qiime", "rescript", "evaluate-fit-classifier",
                "--i-sequences", train_seq,
                "--i-taxonomy", train_tax,
                "--p-n-jobs", threads,
                "--o-classifier", rescript_classifier,
                "--o-evaluation", rescript_evaluation,
                "--o-observed-taxonomy", rescript_observed_taxonomy
            ],
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        logger.error("training failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info("DONE training nb classifier")

    nb_model = outdir / "rescript_classifier.qza"
    nb_classification = outdir / "nb_classification.qza"

    logger.info("running nb classifier")
    try:
        subprocess.run(
            [
                "qiime", "feature-classifier", "classify-sklearn",
                "--i-reads", asv_seqs,
                "--i-classifier", nb_model,
                "--p-n-jobs", threads,
                "--o-classification", nb_classification
            ],
            capture_output=True,
            text=True,
            check=True
        )
    
    except subprocess.CalledProcessError as e:
        logger.error("classification failed")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info("DONE running nb classifier")

    return nb_classification

def run_blast(logger, asv_seqs, ref_tax, ref_seq, params, outdir):
    """Run BLAST on any seqs unclassified after VSEARCH and Naive Bayes."""
    out_tax = outdir / "blast_tax.qza"
    top_hits = outdir / "blast_top_hits.qza"

    logger.info("running blast")
    try:
        subprocess.run(
            [
                "qiime", "feature-classifier", "classify-consensus-blast",
                "--i-query", asv_seqs,
                "--i-reference-reads", ref_seq,
                "--i-reference-taxonomy", ref_tax,
                "--p-perc-identity", params["perc_identity"],
                "--p-query-cov", params["query_cov"],
                "--p-maxaccepts", params["max_accepts"],
                "--o-classification", out_tax,
                "--o-search-results", top_hits
            ],
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        qiime_tmp_dir = e.stderr[-35:-1]
        logger.error(f"blast failed. check {qiime_tmp_dir} for more info.")
    logger.info("DONE running blast")

    return out_tax

def stitch_tax_files(logger, tax_files, outdir):
    """Assemble all taxa retained across all three classification steps and combine them into one tsv."""
    logger.info("assembling final tax file")

    all_tax_dfs = [] # turn each tsv into a dataframe and store it here
    for file in tax_files:
        if file == None:
            logger.error("one file was None type, skipping file")
        current_tax_df = pd.read_csv(file, sep="\t") # read in current taxonomy file as dataframe
        all_tax_dfs.append(current_tax_df)
    final_tax_df = pd.concat(all_tax_dfs).drop_duplicates() # stitch taxonomy files together and drop duplicates
    
    final_tax_tsv = outdir / "final_tax.tsv"
    final_tax_df.to_csv(final_tax_tsv, sep="\t", index=False)

    logger.info(f"DONE assembling tax file")
    return final_tax_tsv

def make_tax_map(logger, tax_file):
    """Read taxonomy file. Store ids with assignments in dict."""
    logger.info("assembling hash for ids and taxonomy")
    with open(tax_file) as f: # read taxonomy file and store each assignment
        f.readline()
        tax_lines = f.readlines()

    tax_hash = {} # map ids to best assignment
    for line in tax_lines:
        otuid, taxon, cons, conf = line.split("\t") # get all elements from the line
        tax_split = taxon.split(";")
        best_assignment = tax_split[-1] # check the level of assignment
        assignment_level = best_assignment[:3] # prefix tells the level of the assignment
        if assignment_level == "s__": # qiime lists both genus and species at the species level
            genus, species = best_assignment[3:].split("_")
            best_assignment = genus + " " + species
        else:
            best_assignment = best_assignment[3:] # get just the name not the level prefix
        tax_hash[otuid] = best_assignment
    
    logger.info("DONE assembling hash")
    return tax_hash

def collapse_unique_hits(logger, feat_tab):
    """Sum reads across duplicates."""
    logger.info("collapsing duplicates and summing counts")
    feat_tab_transpose = feat_tab.T # sets tax assignments to columns
    unique_hits = feat_tab_transpose.columns.unique() # get all hits without duplicates
    unique_cols = {}
    for hit in unique_hits:
        hit_cols = feat_tab_transpose[hit]
        one_hit = len(hit_cols.shape) < 2 # check if there are multiple cols and thus duplicates
        if one_hit:
            unique_cols[hit] = hit_cols # if there's no dupes, just port the data over
        else:
            unique_cols[hit] = hit_cols.sum(axis=1) # sum counts if there's dupes
    
    feat_tab_unique = pd.DataFrame(data=unique_cols).T
    feat_tab_unique.index.rename("Taxon", inplace=True)
    logger.info("DONE getting unique hits")
    return feat_tab_unique

def map_tax_to_feat_table(logger, feat_table, final_tax_tsv, outdir):
    """Map taxonomy assignments to the DADA2 feature table."""
    logger.info("starting tax mapping")
    new_name = "feat_table.biom" # qiime stores dada2 feature table as a biom file
    feat_table_biom = extract_qza(logger, feat_table, new_name, outdir)
    
    logger.info(f"converting {feat_table_biom} to tsv")
    feat_table_tsv = outdir / "feat_table.tsv"
    try:
        subprocess.run( # convert biom file to tsv for parsing
            [
                "biom", "convert",
                "-i", feat_table_biom,
                "-o", feat_table_tsv,
                "--to-tsv"
            ],
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"biom convert failed on {feat_table_biom}")
        logger.error(e.stderr)
        sys.exit(1)
    logger.info(f"converted {feat_table_biom} to {feat_table_tsv}")

    feat_tab_df = pd.read_csv(
        feat_table_tsv,
        sep="\t",
        skiprows=1 # skip the commented row of the converted biom file
    ).rename(columns={"#OTU ID": "Taxon"}) # change name of taxon column
    tax_tab_df = pd.read_csv(final_tax_tsv, sep="\t")

    tax_hash = make_tax_map(logger, final_tax_tsv)

    logger.info("replacing ids with best assignment")
    assigned_tax_col = [] # column to replace ids in final df
    for otuid in feat_tab_df["Taxon"]:
        try:
            assigned_tax_col.append(tax_hash[otuid])
        except KeyError:
            assigned_tax_col.append("Unassigned") # if id was not assigned, mark so that it can be removed later
    
    feat_tab_df["Taxon"] = assigned_tax_col # set taxonomy in place of ids
    feat_tab_assigned = feat_tab_df.set_index("Taxon") # set assignment to index
    feat_tab_assigned_out_tsv = outdir / "feat_tab_assigned.tsv"
    feat_tab_assigned.to_csv(feat_tab_assigned_out_tsv, sep="\t") # write out assigned tsv
    logger.info(f"wrote feat table with best assignment to {feat_tab_assigned_out_tsv}")

    feat_tab_no_ids = feat_tab_assigned.drop("Unassigned", axis=0) # drop rows where the ASV was not assigned
    feat_tab_no_ids_out_tsv = outdir / "feat_tab_no_ids.tsv"
    feat_tab_no_ids.to_csv(feat_tab_no_ids_out_tsv, sep="\t")
    logger.info(f"wrote feat table without unassigned ASVS to {feat_tab_no_ids_out_tsv}")

    feat_tab_unique = collapse_unique_hits(logger, feat_tab_no_ids)
    feat_tab_unique_out_tsv = outdir / "feat_tab_unique.tsv"
    feat_tab_unique.to_csv(feat_tab_unique_out_tsv, sep="\t")
    logger.info(f"wrote feat table without summed dupes to {feat_tab_unique_out_tsv}")

def main():
    config = load_config("config.ini")
    args = get_args(config)
    outdir = create_outdir("results")
    logger = setup_logger(outdir / "pipeline.log")

    threads = args.threads
    logger.info(f"starting pipeline with {threads} threads")

    if not args.skip_to_mapping:
        crabs_ref_tax = Path(args.crabs_database_tax).resolve()
        crabs_ref_seq = Path(args.crabs_database_seq).resolve()
        if not crabs_ref_tax.is_file():
            logger.error(f"database file {crabs_ref_tax} not found")
            sys.exit(1)
        if not crabs_ref_seq.is_file():
            logger.error(f"database file {crabs_ref_seq} not found")
            sys.exit(1)
        logger.info("loaded crabs reference db files")

        blast_ref_tax = Path(args.blast_database_tax).resolve()
        blast_ref_seq = Path(args.blast_database_seq).resolve()
        if not blast_ref_tax.is_file():
            logger.error(f"database file {blast_ref_tax} not found")
            sys.exit(1)
        if not blast_ref_seq.is_file():
            logger.error(f"database file {blast_ref_seq} not found")
            sys.exit(1)
        logger.info("loaded blast reference db files")

        if not args.archive: # only import reads if archive not supplied (this takes a while)
            reads_archive = Path(outdir / "reads.qza").resolve()
            if not args.manifest:
                print("please supply manifest file using the --manifest option")
                logger.error("no manifest file supplied")
                sys.exit(1)
            manifest = Path(args.manifest).resolve() # only load manifest if there's no archive supplied
            if not manifest.is_file():
                logger.error(f"cannot access manifest {manifest}. file does not exist.")
                sys.exit(1)
            logger.info(f"loaded manifest: {manifest}")
            import_reads(logger, manifest, reads_archive)
        else:
            logger.info("skipping import and using archive file instead")
            reads_archive = Path(args.archive).resolve() # get archive path
        logger.info(f"loaded archive: {reads_archive}")

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

        vsearch_dir = outdir / "vsearch"
        vsearch_dir.mkdir()
        vsearch_out = run_vsearch(logger, asv_seqs, crabs_ref_tax, crabs_ref_seq, threads, vsearch_dir)

        species_level = 7 # looking for exact matches all the way to species level
        vsearch_unassigned_tax, vsearch_retained_tax = parse_output(logger, vsearch_out, vsearch_dir, species_level)
        vsearch_unassigned_seq_archive, vsearch_unassigned_seq_fasta = filter_seqs(logger, asv_seqs, vsearch_unassigned_tax, vsearch_dir)

        bayes_dir = outdir / "bayes"
        bayes_dir.mkdir()
        bayes_out = run_nb_classifier(logger, vsearch_unassigned_seq_archive, crabs_ref_tax, crabs_ref_seq, threads, bayes_dir)

        family_level = 5 # looking only above family level
        bayes_unassigned_tax, bayes_retained_tax = parse_output(logger, bayes_out, bayes_dir, family_level)
        bayes_unassigned_seq_archive, bayes_unassigned_seq_fasta = filter_seqs(logger, asv_seqs, bayes_unassigned_tax, bayes_dir)

        blast_params = { # blast parameters from command line arguments
            "perc_identity": args.perc_identity,
            "query_cov": args.query_cov,
            "max_accepts": args.max_accepts,
            "threads": threads
        }
        blast_dir = outdir / "blast"
        blast_dir.mkdir()
        blast_out = run_blast(logger, bayes_unassigned_seq_archive, blast_ref_tax, blast_ref_seq, blast_params, blast_dir)

        if blast_out.is_file(): # only run these if blast completes successfully
            blast_unassigned_tax, blast_retained_tax = parse_output(logger, blast_out, blast_dir, family_level)
            blast_unassigned_seq_archive, blast_unassigned_seq_fasta = filter_seqs(logger, asv_seqs, blast_unassigned_tax, bayes_dir)
        else:
            blast_retained_taxa = None
            blast_unassigned_seq_fasta = None

        tax_files = [
            vsearch_retained_tax,
            bayes_retained_tax,
            blast_retained_tax
        ]
    
    else:
        logger.info("skipping straight to mapping step")
        tax_files = [
            "post-tax-files/vsearch_retained_tax.tsv",
            "post-tax-files/nb_retained_tax.tsv",
            "post-tax-files/blast_retained_tax.tsv"
        ]
        feat_table = "post-tax-files/feat_table.qza"

    mapping_dir = outdir / "mapping_files"
    mapping_dir.mkdir()
    final_tax_tsv = stitch_tax_files(logger, tax_files, mapping_dir)
    map_tax_to_feat_table(logger, feat_table, final_tax_tsv, mapping_dir)

    logger.info("pipeline end")

if __name__ == "__main__":
    main()