# use metadata to find paths to read file for each sample and generate manifest for each replicate

import sys
import os
import pandas as pd
import argparse
from pathlib import Path
import re

def parse_args():
    """
    parse arguments, making sure data is accessible and study name is correct
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metadata", required=True, help="metadata file mapping sample ids to study")
    parser.add_argument("-s", "--study", required=True, help="possible study types: DamBaseline, JuneJulyTemporal, EbonyTemporal, Filter")
    parser.add_argument("-d", "--data", required=True, help="directory with read files to be searched recursively")
    parser.add_argument("-o", "--manifest", help="name of generated manifest file")
    args = parser.parse_args()

    studies = ["DamBaseline", "JuneJulyTemporal", "EbonyTemporal", "Filter"]
    if args.study not in studies:
        print("possible studies include:", end=" ")
        for study in studies:
            print(study, end=" ")
        print()
        sys.exit(1)

    if not os.path.exists(args.metadata): # make sure metadata is accessible
        print("error: metadata not accessible")
        sys.exit(1)

    if not os.path.isdir(args.data): # make sure data is accessible
        print("error: disk not mounted")
        sys.exit(1)
    
    return args

def subset_data(metadata, study):
    """
    subset data for samples from target study
    """
    meta_df = pd.read_excel(metadata) # read in metadata
    meta_subset = meta_df[~meta_df[study].isna()] # subset for just target study
    subset_sams = list(meta_subset["Sample ID"]) # grab sample ids
    with open(f"{study}_ids.txt", "w") as f:
        for sid in subset_sams:
            f.write(f"{sid}\n")
    return subset_sams

def generate_manifest(data, study, manifest, sam_ids):
    """
    grab read data from target study and write paths out to manifest tsv
    """
    if not manifest:
        manifest = f"{study}_manifest.tsv"
    if not manifest.endswith(".tsv"):
        manifest = manifest + ".tsv"

    exclude_dirs = {"trimmed", "mussel", "$RECYCLE.BIN", "extra", "Picq04_4.15.2026"} # ignore these directories for now
    read_paths = [ # get all fastq files from all subdirs
        str(p) for p in Path(data).rglob("*.fastq.gz") # recursively extract all fastq files
        if exclude_dirs.isdisjoint(p.parts) # exclude paths that include these directories
    ]

    seen_ids = [] # store ids that have already been written to catch dupes
    with open(manifest, "w") as f:
        f.write(f"sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for i in range(0, len(read_paths), 2): # grab every other file since the reads are paired
            forward_path = read_paths[i]
            for sid in sam_ids:
                pattern = rf"[/_]{sid}D?[-_]" # regex to match. accounts for ids preceeded by SP-## and ids that contain other ids
                if re.search(pattern, forward_path): # only grab file if it's from the target study
                    print(sid, forward_path)
                    reverse_path = read_paths[i+1]
                    split_file_name = os.path.split(read_paths[i])[1].split("_")
                    rep_id = split_file_name[1] if split_file_name[0].startswith("SP") else split_file_name[0] # get sample id with rep number
                    if rep_id.endswith("mussel"): # skip mussel files
                        continue
                    f.write(f"{rep_id}\t{forward_path}\t{reverse_path}\n")
                seen_ids.append(sid)
    return

def amend_metadata():
    """
    add 
    """
    pass

def main():
    args = parse_args()
    sam_ids = subset_data(args.metadata, args.study)
    generate_manifest(args.data, args.study, args.manifest, sam_ids)

if __name__ == "__main__":
    main()