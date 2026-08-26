# PRE-PROCESSING: find read files for target study, generate manifest and metadata for each replicate

import sys
import os
import pandas as pd
import argparse
from pathlib import Path
import re

def parse_args():
    """parse arguments, making sure data is accessible and study name is correct"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metadata", help="metadata file mapping sample ids to study")
    parser.add_argument("-s", "--study", help="possible study types: DamBaseline, JuneJulyTemporal, EbonyTemporal, Filter_5.0v0.45")
    parser.add_argument("-d", "--data", required=True, help="directory with read files to be searched recursively")
    parser.add_argument("-b", "--blanks", action="store_true", help="generate manifest for extraction and field blanks and positive and negative samples only")
    args = parser.parse_args()

    if not os.path.isdir(args.data): # make sure data is accessible
        print("error: disk not mounted")
        sys.exit(1)
    
    return args

def subset_data(metadata, study):
    """subset data for samples from target study"""
    meta_df = pd.read_excel(metadata) # read in metadata
    
    blank_cols = meta_df.columns[meta_df.columns.str.startswith("Unnamed")] # find metadata cols that were blank in excel
    meta_df = meta_df.drop(blank_cols, axis=1) # drop those cols
    meta_df = meta_df[~meta_df[study].isna()] # subset for study

    subset_sams = list(meta_df["Sample ID"]) # grab sample ids
    meta_df.index = meta_df["Sample ID"] # set rownames to sample id

    return meta_df, subset_sams

def generate_manifest(data, study, sam_ids):
    """grab read data from target study and write paths out to manifest tsv"""
    manifest = f"{study}_manifest.tsv"

    exclude_dirs = {"trimmed", "mussel", "$RECYCLE.BIN", "extra"} # ignore these directories for now
    read_paths = [ # get all fastq files from all subdirs
        str(p) for p in Path(data).rglob("*.fastq.gz") # recursively extract all fastq files
        if exclude_dirs.isdisjoint(p.parts) # exclude paths that include these directories
    ]

    manif_ids = [] # store ids used to generate manifest
    with open(manifest, "w") as f:
        f.write(f"sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for i in range(0, len(read_paths), 2): # grab every other file since the reads are paired
            forward_path = read_paths[i]
            for sid in sam_ids:
                pattern = rf"[/_]{sid}D?[-_]" # regex to match. accounts for ids preceeded by SP-## and ids that contain other ids
                if re.search(pattern, forward_path): # only grab file if it's from the target study
                    reverse_path = read_paths[i+1]
                    split_file_name = os.path.split(read_paths[i])[1].split("_")
                    rep_id = split_file_name[1] if split_file_name[0].startswith("SP") else split_file_name[0] # get sample id with rep number
                    if rep_id.endswith("mussel"): # skip mussel files
                        continue
                    if rep_id in manif_ids: # if we've already seen this sample, skip it (helps with metadata dupes)
                        continue
                    f.write(f"{rep_id}\t{forward_path}\t{reverse_path}\n")
                    manif_ids.append(rep_id)

    return manif_ids

def match_ids(final_df, meta_df, rep_id):
    """resolve duplcates and mismatched ids"""
    meta_sample_id = "-".join(rep_id.split("-")[:2])
    in_metadata = meta_sample_id in meta_df.columns # if the sample id is in the metadata (not a D, i.e.)
    if not in_metadata: # if the input sample id is not in the metadata
        if meta_sample_id.endswith("D"): # if the sample id is a duplicate
            meta_sample_id = meta_sample_id[:-1] # index for the non-duplicate id
        else:
            matched_line = ["missing"]*len(meta_df) # if its simply not in the metadata, note that
            return final_df # and return
    
    duplicate = len(meta_df[meta_sample_id].shape) > 1 # if there is more than one column, its a duplicate
    if duplicate:
        matched_line = meta_df[meta_sample_id].iloc[:, 0] # grab first col only
    else: # if there are no problems
        matched_line = meta_df[meta_sample_id] # just port data over

    return matched_line

def generate_metadata(manif_ids, meta_df, study):
    """add sample metadata to replicates"""
    t_meta_df = meta_df.T # transpose metadata
    final_merged_df = pd.DataFrame() # initialize df to output
    final_merged_df.index = t_meta_df.index

    final_cols = {}
    for rep_id in manif_ids:
        final_cols[rep_id] = match_ids(final_merged_df, t_meta_df, rep_id)

    final_df = pd.DataFrame(final_cols).T
    final_df.index = final_df.index.rename("Replicate ID")
    final_df.to_csv(f"{study}_metadata.tsv", sep="\t") # un-transpose and write out

    return

def get_blanks(data):
    """Generate manifest file for just extraction blanks, etc."""
    manifest = f"blank_manifest.tsv"

    pattern = re.compile(r"(^EB-|FB-|^pos|^neg)") # match this regex to get extraction blanks
    exclude_dirs = {"trimmed", "mussel", "$RECYCLE.BIN", "extra", "Picq04_4.15.2026"} # ignore these directories for now
    read_paths = [ # get all fastq files from all subdirs
        str(p) for p in Path(data).rglob("*.fastq.gz") # recursively extract all fastq files
        if exclude_dirs.isdisjoint(p.parts) # exclude paths that include these directories
        if pattern.search(p.name)
    ]

    sample_ids = []
    with open(manifest, "w") as f:
        f.write(f"sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for i in range(0, len(read_paths), 2): # grab every other file since the reads are paired
            forward_path = read_paths[i]
            reverse_path = read_paths[i+1]
            split_file_name = os.path.split(read_paths[i])[1].split("_")
            rep_id = split_file_name[1] if split_file_name[0].startswith("SP") else split_file_name[0] # get sample id with rep number
            if rep_id in sample_ids: # don't add duplicates
                continue
            if rep_id.endswith("mussel"): # ignore mussel blanks for now
                continue
            f.write(f"{rep_id}\t{forward_path}\t{reverse_path}\n")
            sample_ids.append(rep_id)

    return

def main():
    args = parse_args()
    if not args.blanks:
        if not args.metadata:
            print("Provide metadata file")
            sys.exit(1)
        if not args.study:
            print("Provide study")
            sys.exit(1)
        meta_df, sam_ids = subset_data(args.metadata, args.study)
        manif_ids = generate_manifest(args.data, args.study, sam_ids)
        generate_metadata(manif_ids, meta_df, args.study)
    else:
        get_blanks(args.data)

if __name__ == "__main__":
    main()