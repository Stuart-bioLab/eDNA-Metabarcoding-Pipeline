# Collect samples. Generate manifest and blank map.abs

import sys
import argparse
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict
import shutil
import numpy as np

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", help="Path to dir containing sequence read files.")
    parser.add_argument("-s", "--study", help="Target study to subset reads for.")
    parser.add_argument("-f", "--file", help="Provide tsv mapping read prefix to filepath.")
    args = parser.parse_args()
    return args

def subset_metadata(metadata, study, outdir):
    """Subset metadata for target study. Return subset DataFrame and sample names for study."""
    meta_df = pd.read_excel(metadata) # read in excel file
    unnamed_cols = meta_df.columns[meta_df.columns.str.startswith("Unnamed")] # find unnamed columns (i don't know why they're here)
    meta_df.drop(unnamed_cols, axis=1, inplace=True) # drop unnamed cols

    meta_subset_df = meta_df[~meta_df[study].isna()] # subset dataframe for samples in target study
    meta_subset_df.to_csv(outdir / f"{study}_subset_metadata.tsv", sep="\t", index=False) # write out maybe i'll want to look at this idk

    sam_ids = list(meta_subset_df["Sample ID"]) # get list of sample names from study
    no_nan_sam_ids = [x for x in sam_ids if x is not np.nan]

    return meta_df, meta_subset_df, no_nan_sam_ids

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

def group_by_field_blank(df, meta_df, study, outdir):
    """Write out file showing which field blank each sample was collected with."""
    df_groupby_date = df.groupby("Date Collected")["Sample ID"].apply(list) # attribute a list of samples collected to each date
    date_dict = df_groupby_date.to_dict() # convert to a dict where timestamp: [samples collected]

    outfile = outdir / f"{study}_field_blank_map.tsv"
    with open(outfile, "w") as f:
        f.write(f"sample-id\t{study}-field-blank-id\tdate-collected\n")
        written_ids = []
        for k, v in date_dict.items():
            no_nan_sams = [x for x in v if x is not np.nan]
            fb = [x for x in no_nan_sams if "FB" in x]
            if not fb: # if field blank is not in DamBaseline subset
                all_samples_that_day = meta_df[meta_df["Date Collected"] == k]["Sample ID"].dropna() # look at all samples collected on target date
                fb = [x for x in all_samples_that_day if "FB" in x] # isolate the field blank
            date_collected = f"{k.month}/{k.day}/{k.year}" # reformat date
            for sam_id in no_nan_sams:
                if len(fb) > 1:
                    for fblank in fb:
                        if sam_id[-2:] == fblank[-2:]:
                            single_field_blank = fblank
                else:
                    single_field_blank = "".join(fb) if fb else "Not sequenced"
                if sam_id is not single_field_blank:
                    if sam_id not in written_ids:
                        f.write(f"{sam_id}\t{single_field_blank}\t{date_collected}\n")
                        written_ids.append(sam_id)

    return outfile

def get_filepaths(study_ids, infile):
    """Create dict that groups all read filepaths with their corresponding fastq prefix."""
    with open(infile, "r") as f:
        f.readline()
        lines = f.readlines()

    filepath_dict = defaultdict(list) # prefix: [read filepaths]
    for l in lines:
        prefix, filepath = l.split("\t")
        if "mussel" in filepath: # not looking at mussel sams right now
            continue
        if prefix in study_ids:
            filepath_dict[prefix].append(filepath.strip())

    return filepath_dict

def write_sample_manifest(study_ids, reads_list, study, outdir):
    """Generate manifest file containing read filepaths for all available samples in target study"""
    filepath_dict = get_filepaths(study_ids, reads_list)

    outfile = outdir / f"{study}_manifest.tsv"
    seen_rep_ids = [] # store replicate ids here so we don't enter duplicates into the manifest
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
                if replicate_id not in seen_rep_ids:
                    f.write(f"{replicate_id}\t{forward_path}\t{reverse_path}\n")
                    seen_rep_ids.append(replicate_id)
    
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

def build_extraction_map(infile, id_list, study, outdir):
    """Write file mapping each sample to its extraction blank. Then, build map from sample id to eblank id for building manifest."""
    df = pd.read_excel(infile)
    target_samples = df["Sample ID"].isin(id_list) # only get samples from the target study
    subset_df = df[target_samples]
    extr_blank_map_df = subset_df[["Sample ID", "Extraction Negative"]] # just want sample ID and its eblank

    outfile = outdir / f"{study}_extraction_blank_map.tsv"
    extr_blank_map_df.to_csv(outfile, index=False, sep="\t")

    extr_blank_rep_ids = list( # get a list of extraction blanks
        extr_blank_map_df["Extraction Negative"] # subset for just blank replicate ids
        .dropna() # get rid of blank entries
    )

    eblank_sam_ids = []
    for rep_id in extr_blank_rep_ids:
        sam_id = re.sub("-rep.*$", "", rep_id) # drop replicate number
        if sam_id not in eblank_sam_ids:
            eblank_sam_ids.append(sam_id)

    return eblank_sam_ids, outfile

def append_extraction_blanks(input_manif, eblank_metadata, reads_list, study, outdir):
    """Append extraction blanks to manifest so they can be run alongside regular samples."""
    rep_id_list = read_replicate_ids(input_manif)
    sam_list, eblank_map = build_extraction_map(eblank_metadata, rep_id_list, study, outdir)

    with open(reads_list, "r") as f:
        f.readline()
        lines = [x.strip() for x in f.readlines()]
    
    final_manif = outdir / f"final_{study}_manifest.tsv"
    shutil.copy(input_manif, final_manif)

    seen_rep_ids = [] # store rep ids here just like with samples
    with open(final_manif, "a") as m:
        for i in range(0, len(lines), 2): # reads are paired
            sam_id, fpath = lines[i].split("\t")
            rpath = lines[i+1].split("\t")[1]
            if sam_id in sam_list: # grab files from target study
                split_path_name = fpath.split("/")[-1].split("_")
                rep_id = split_path_name[1] if split_path_name[0].startswith("SP") else split_path_name[0] # get prefix from filepath (accounting for NWern sams that have SP- in front)
                if "mussel" in rep_id: # not looking a mussel sams right now
                    continue
                fpath = fpath.replace("/mnt/d/", "/mnt/g/") # data on the PC is stored in G:
                rpath = rpath.replace("/mnt/d/", "/mnt/g/")
                if rep_id not in seen_rep_ids: # don't repeat replicates
                    m.write(f"{rep_id}\t{fpath}\t{rpath}\n")
                    seen_rep_ids.append(rep_id)

    return eblank_map

def build_replicate_metadata():
    pass

def main():
    args = get_args()
    data = args.data
    reads_list = args.file
    study = args.study
    metadata = "LIVE_FoxRiver_eDNA_Field_Data_Clean.xlsx"
    blank_map = "all_sample_metadata.xlsx"

    studies = ["DamBaseline", "JuneJulyTemporal", "EbonyTemporal", "Filter_5.0v0.45"]

    if study == "d":
        study = "DamBaseline"
    if study == "j":
        study = "JuneJulyTemporal"
    if study == "e":
        study = "EbonyTemporal"
    if study == "f":
        study = "Filter_5.0v0.45"

    if study not in studies:
        print("Possible studies include", end=" ")
        print(", ".join(studies))
        sys.exit(1)

    outdir = Path("sample_mapping_files") / study
    outdir.mkdir(exist_ok=True, parents=True)
    full_meta_df, meta_study_df, study_ids = subset_metadata(metadata, study, outdir)
    
    if data:
        reads_list = find_read_files(data, outdir)
    else:
        if not reads_list:
            print("Supply path to data dir.")
            sys.exit(1)
    
    fblank_map_file = group_by_field_blank(meta_study_df, full_meta_df, study, outdir)

    study_manif = write_sample_manifest(study_ids, reads_list, study, outdir)
    eblank_map_file = append_extraction_blanks(study_manif, blank_map, reads_list, study, outdir)

    build_replicate_metadata(fblank_map_file, eblank_map_file)

if __name__ == "__main__":
    main()