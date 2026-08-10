# append sample metadata to replicates given metadata and manifest files

import pandas as pd
import argparse

def get_args():
    """
    get arguments from command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--metadata", required=True, help="metadata file mapping sample ids to study")
    parser.add_argument("-s", "--study", required=True, help="possible study types: DamBaseline, JuneJulyTemporal, EbonyTemporal, Filter")
    parser.add_argument("-o", "--manifest", required=True, help="name of generated manifest file")
    return parser.parse_args()

def import_data(meta, mani, study):
    """
    import and process data
    """
    meta_df = pd.read_excel(meta)
    blank_cols = meta_df.columns[meta_df.columns.str.startswith("Unnamed")] # find metadata cols that were blank in excel
    meta_df = meta_df.drop(blank_cols, axis=1) # drop those cols
    meta_df = meta_df[~meta_df[study].isna()] # subset for study
    meta_df.index = meta_df["Sample ID"] # set rownames to sample id
    meta_df = meta_df.drop("Sample ID", axis=1) # and drop

    mani_df = pd.read_csv(mani, sep="\t")
    mani_df["metadata_sample_id"] = ["-".join(x.split("-")[:2]) for x in mani_df["sample-id"]] # add column without replicate number for mapping to metadata

    return meta_df, mani_df

def match_ids(final_df, meta_df, rep_id):
    """
    resolve duplcates and mismatched ids
    """
    meta_sample_id = "-".join(rep_id.split("-")[:2])
    in_metadata = meta_sample_id in meta_df.columns # if the sample id is in the metadata (not a D, i.e.)
    if not in_metadata: # if the input sample id is not in the metadata
        if meta_sample_id.endswith("D"): # if the sample id is a duplicate
            meta_sample_id = meta_sample_id[:-1] # index for the non-duplicate id
        else:
            final_df[rep_id] = ["missing"]*len(meta_df) # if its simply not in the metadata, note that
            return final_df # and return
    
    duplicate = len(meta_df[meta_sample_id].shape) > 1 # if there is more than one column, its a duplicate
    if duplicate:
        final_df[rep_id] = meta_df[meta_sample_id].iloc[:, 0] # grab first col only
    else: # if there are no problems
        final_df[rep_id] = meta_df[meta_sample_id] # just port data over
    return final_df

def merge_dfs(meta_df, mani_df):
    """
    match data from metadata to replicate ids from manifest
    """
    t_meta_df = meta_df.T # transpose metadata

    final_merged_df = pd.DataFrame() # initialize df to output
    final_merged_df.index = t_meta_df.index

    for rep_id in mani_df["sample-id"]:
        final_merged_df = match_ids(final_merged_df, t_meta_df, rep_id)

    return final_merged_df.T # un-transpose and return

def main():
    args = get_args()
    meta_df, mani_df = import_data(args.metadata, args.manifest, args.study)
    final_merged_df = merge_dfs(meta_df, mani_df)
    final_merged_df.to_csv("db_merged_meta.tsv", sep="\t")

if __name__ == "__main__":
    main()