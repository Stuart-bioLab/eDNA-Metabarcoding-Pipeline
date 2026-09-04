# Collect samples. Generate manifest and blank map.abs

import sys
import argparse
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict
import shutil

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", help="Path to dir containing sequence read files.")
    parser.add_argument("-s", "--study", help="Target study to subset reads for.")
    parser.add_argument("-f", "--file", help="Provide tsv mapping read prefix to filepath.")
    parser.add_argument("-a", "--add_reads", help="Comma delimited list of prefixes to also add to manifest. Add these after manually finding typos etc.")
    args = parser.parse_args()
    return args

def subset_metadata(metadata, study, outdir):
    """Subset metadata for target study. Return subset DataFrame and sample names for study."""
    meta_df = pd.read_excel(metadata) # read in excel file
    unnamed_cols = meta_df.columns[meta_df.columns.str.startswith("Unnamed")] # find unnamed columns (i don't know why they're here)
    meta_df.drop(unnamed_cols, axis=1, inplace=True) # drop unnamed cols

    meta_subset_df = meta_df[~meta_df[study].isna()] # subset dataframe for samples in target study
    meta_subset_df.to_csv(outdir / "subset_metadata.tsv", sep="\t", index=False) # write out maybe i'll want to look at this idk

    sam_ids = list(meta_subset_df["Sample ID"]) # get list of sample names from study

    return meta_df, meta_subset_df, sam_ids

def find_read_files(data, outdir):
    """Search input dir for fastq files. Write out each read file paired with its prefix."""
    exclude_dirs = {"trimmed", "mussel", "$RECYCLE.BIN", "extra", "Picq04_4.15.2026", "mitogenome_extra", "trimmed_fastq"} # ignore these directories for now
    read_paths = [ # get all fastq files from all subdirs
        p.resolve() for p in sorted(Path(data).rglob("*.fastq.gz")) # recursively extract all fastq files
        if exclude_dirs.isdisjoint(p.parts) # exclude paths that include above directories
    ]

    outfile = outdir / "read_filepath_list.tsv"
    with open(outfile, "w") as f:
        f.write(f"read-prefix\tread-filepath\n") # header
        for p in read_paths: # iterate over all fastq files found by rglob
            split_path = p.name.split("_") # split just the name of the read file
            prefix = split_path[1] if split_path[0].startswith("SP") else split_path[0] # get prefix, dropping SP-## if necessary
            prefix_no_rep = re.sub("-rep.*", "", prefix) # remove replicate number from prefix so its easier to match with metadata
            f.write(f"{prefix_no_rep}\t{str(p)}\n") # write prefix and filepath
    
    return outfile

def parse_reads_list(infile):
    """Read in all unique read path prefixes."""
    with open(infile, "r") as f:
        f.readline()
        prefix_list = []
        for line in f.readlines():
            sam_prefix = line.split("\t")[0]
            if sam_prefix not in prefix_list:
                prefix_list.append(sam_prefix)
    
    return prefix_list

def write_out_dict(sam_dict, dupes_list, outdir):
    """Write out map to facilitate manual disambiguation. Also, used next to grab all possible read files for manifest."""
    dict_out = outdir / "sample_id_best_matches.tsv"
    with open(dict_out, "w") as f:
        f.write("sample-id\tpossible-read-files\tno-read-file-found\n")
        for k, v in sam_dict.items():
            candidate_read_files = ";".join(v)
            if candidate_read_files == "NA":
                f.write(f"{k}\t\t{candidate_read_files}\n")
            else:
                f.write(f"{k}\t{candidate_read_files}\n")

    dupes_out = outdir / "sample_dupes.txt"
    with open(dupes_out, "w") as f:
        for d in dupes_list:
            f.write(f"{d}\n")

    return dict_out

def match_sample_ids(infile, meta_sam_ids, outdir):
    """Match sequence read prefix with sample names. If samples do not have the same name, find the next best match."""
    read_prefixes = parse_reads_list(infile)
    
    sample_id_map = defaultdict(list)
    duplicates = [] # store dupes here
    for s in meta_sam_ids: # iterate over samples in the metadata
        if s not in sample_id_map.keys(): # i.e., if this is not a duplicate
            if s in read_prefixes: # if there is a 1-to-1 map from metadata to read files
                sample_id_map[s].append(s) # the names are identical
            elif s[:-1] in read_prefixes: # check if samples were pooled, e.g. NAD-DNL1 and NAD-DNL2 pooled as NAD-DNL
                sample_id_map[s].append(s[:-1])
            else: # otherwise, find close file names
                pattern = rf"{s[:-1]}.*" # check for typos
                for p in read_prefixes:
                    if re.match(pattern, p):
                        sample_id_map[s].append(p)
            for p in read_prefixes: # also check for possible duplicates
                pattern = rf"{s}D" # duplicates have the same file name but with a D appended
                if re.match(pattern, p):
                    sample_id_map[s].append(p)
        else:
            duplicates.append(s)
        if s not in sample_id_map.keys():
            sample_id_map[s].append("NA")

    sample_read_map = write_out_dict(sample_id_map, duplicates, outdir)
    return sample_read_map    

def group_by_field_blank(df, meta_df, study, outdir):
    """Write out file showing which field blank each sample was collected with."""
    df_groupby_date = df.groupby("Date Collected")["Sample ID"].apply(list) # attribute a list of samples collected to each date
    date_dict = df_groupby_date.to_dict() # convert to a dict where timestamp: [samples collected]

    outfile = outdir / "field_blank_map.tsv"
    with open(outfile, "w") as f:
        f.write(f"sample-id\t{study}-field-blank\tother-field-blank\tdate-collected\n")
        non_study_fbs = []
        for k, v in date_dict.items():
            study_fb = None # if field blank is on current study
            other_fb = None # if field blank is on a different study
            study_fb = [x for x in v if "FB" in x]
            if not study_fb: # if field blank is not in target study subset (was collected on same day with other study)
                all_samples_that_day = meta_df[meta_df["Date Collected"] == k]["Sample ID"] # look at all samples collected on target date
                other_fb = [x for x in all_samples_that_day if "FB" in x] # isolate the field blank
                for fb in other_fb: # add field blank even if its not on study
                    non_study_fbs.append(fb)
            date_collected = f"{k.month}/{k.day}/{k.year}" # reformat date
            study_fb = "NA" if not study_fb else "".join(study_fb)
            other_fb = "NA" if not other_fb else "".join(other_fb)
            for sam_id in v:
                f.write(f"{sam_id}\t{study_fb}\t{other_fb}\t{date_collected}\n")
    
    if not non_study_fbs:
        non_study_fbs = None

    return outfile, non_study_fbs

def get_read_file_prefixes(infile):
    """Get all unique filepath prefixes that correspond to metadata entries."""
    read_prefix_list = []
    with open(infile, "r") as f:
        f.readline()
        for line in f.readlines():
            read_files = line.split("\t")[1].split(";")
            for r in read_files:
                r = r.strip()
                if not r: # if there's no entry
                    continue
                if r not in read_prefix_list:
                    read_prefix_list.append(r)

    return read_prefix_list

def get_filepaths(prefix_list, in_file):
    """Create dict that groups all read filepaths with their corresponding fastq prefix."""
    with open(in_file, "r") as f:
        f.readline()
        lines = f.readlines()

    filepath_dict = defaultdict(list) # prefix: [read filepaths]
    for l in lines:
        prefix, filepath = l.split("\t")
        if "mussel" in filepath: # not looking at mussel sams right now
            continue
        if prefix in prefix_list:
            filepath_dict[prefix].append(filepath.strip())

    return filepath_dict

def write_sample_manifest(sample_read_map, reads_list, study, other_reads, outdir):
    """Generate manifest file containing read filepaths for all available samples in target study"""
    read_prefix_list = get_read_file_prefixes(sample_read_map)
    if other_reads:
        for r in other_reads: # add manually input reads
            if r not in read_prefix_list:
                read_prefix_list.append(r)
    filepath_dict = get_filepaths(read_prefix_list, reads_list)

    outfile = outdir / f"{study}_manifest.tsv"
    with open(outfile, "w") as f:
        f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for v in filepath_dict.values():
            sorted_filepaths = []
            for filepath in sorted(v):
                sorted_filepaths.append(filepath.replace("/mnt/d/", "/mnt/g/")) # data is mounted on a different drive on the PC
            for i in range(0, len(sorted_filepaths), 2):
                forward_path = sorted_filepaths[i]
                reverse_path = sorted_filepaths[i+1]
                split_file_name = forward_path.split("/")[-1].split("_")
                replicate_id = split_file_name[1] if split_file_name[0].startswith("SP") else split_file_name[0]
                f.write(f"{replicate_id}\t{forward_path}\t{reverse_path}\n")
    
    return outfile

def read_replicate_ids(infile):
    """Extract all replicate IDs from sample manifest."""
    replicate_id_list = []
    with open(infile, "r") as f:
        f.readline()
        for line in f.readlines():
            repid = line.split("\t")[0]
            if repid not in replicate_id_list:
                replicate_id_list.append(repid)

    return replicate_id_list

def build_extraction_map(infile, id_list, outdir):
    """Write file mapping each sample to its extraction blank."""
    df = pd.read_excel(infile)
    target_samples = df["Sample ID"].isin(id_list) # only get samples from the target study
    subset_df = df[target_samples]
    extr_blank_map_df = subset_df[["Sample ID", "Extraction Negative"]] # just want sample ID and its eblank

    outfile = outdir / "extraction_blank_map.tsv"
    extr_blank_map_df.to_csv(outfile, index=False, sep="\t")

    extr_blank_rep_ids = list( # get a list of extraction blanks
        extr_blank_map_df["Extraction Negative"] # subset for just blank replicate ids
        .dropna() # get rid of blank entries
    )

    eblank_sam_ids = []
    for rep_id in extr_blank_rep_ids:
        sam_id = re.sub("-rep.*$", "", rep_id)
        if sam_id not in eblank_sam_ids:
            eblank_sam_ids.append(sam_id)

    return eblank_sam_ids

def append_extraction_blanks(input_manif, eblank_metadata, reads_list, study, outdir):
    """Append extraction blanks to manifest so they can be run alongside regular samples."""
    rep_id_list = read_replicate_ids(input_manif)
    sam_list = build_extraction_map(eblank_metadata, rep_id_list, outdir)

    with open(reads_list, "r") as f:
        f.readline()
        lines = [x.strip() for x in f.readlines()]
    
    final_manif = outdir / f"final_{study}_manifest.tsv"
    shutil.copy(input_manif, final_manif)

    with open(final_manif, "a") as m:
        for i in range(0, len(lines), 2):
            sam_id, fpath = lines[i].split("\t")
            rpath = lines[i+1].split("\t")[1]
            if sam_id in sam_list:
                split_path_name = fpath.split("/")[-1].split("_")
                rep_id = split_path_name[1] if split_path_name[0].startswith("SP") else split_path_name[0]
                if "mussel" in rep_id:
                    continue
                fpath = fpath.replace("/mnt/d/", "/mnt/g/")
                rpath = rpath.replace("/mnt/d/", "/mnt/g/")
                m.write(f"{rep_id}\t{fpath}\t{rpath}\n")

def main():
    args = get_args()
    data = args.data
    reads_list = args.file
    study = args.study
    metadata = "FReDNA_master_metadata.xlsx"
    blank_map = "all_sample_metadata.xlsx"
    input_reads = args.add_reads

    studies = ["DamBaseline", "JuneJulyTemporal", "EbonyTemporal", "Filter_5.0v0.45"]
    if study == "Filter":
        study = "Filter_5.0v0.45"
    if study not in studies:
        print("Possible studies include", end=" ")
        print(", ".join(studies))
        sys.exit(1)

    outdir = Path("metadata_readfile_bridge") / study
    outdir.mkdir(exist_ok=True, parents=True)
    full_meta_df, meta_study_df, study_ids = subset_metadata(metadata, study, outdir)
    
    if data:
        reads_list = find_read_files(data, outdir)
    else:
        if not reads_list:
            print("Supply path to data dir.")
            sys.exit(1)
    
    sample_read_map = match_sample_ids(reads_list, study_ids, outdir)
    field_blank_map, non_study_fbs = group_by_field_blank(meta_study_df, full_meta_df, study, outdir)

    if non_study_fbs: # add field blanks that are not on target study to manifest so they can also be processed (if files exist)
        if not input_reads: # if there are no other prefixes given at command line
            other_reads = non_study_fbs
        else: # otherwise, search for other field blanks if files exist
            other_reads = input_reads.split(",")
            for fb in non_study_fbs: # set field blanks up to be added to manifest if they are non on target study
                other_reads.append(fb)

    study_manif = write_sample_manifest(sample_read_map, reads_list, study, other_reads, outdir)
    append_extraction_blanks(study_manif, blank_map, reads_list, study, outdir)

if __name__ == "__main__":
    main()