import pandas as pd
from pathlib import Path
from collections import defaultdict
import re
import itertools

def sum_replicates(feat_tab):
    """Split input dataframe into samples with replicates and singletons. Sum read counts across replicates. Concatenate and return dereplicated dataframe."""
    rep_ids = feat_tab.columns
    replicated_sams = rep_ids[rep_ids.str.contains("rep")] # pull out which samples have replicates
    feat_tab_reps = feat_tab[replicated_sams] # get samples with replicates
    feat_tab_non_reps = feat_tab.drop(replicated_sams, axis=1) # get samples without replicates

    sam_ids = defaultdict(list)
    for rep in feat_tab_reps.columns: # map replicates to sample prefixes
        sam = rep[:-5]
        sam_ids[sam].append(rep) # store sample id with list of its replicates

    final_tab_cols = {}
    for sam, reps in sam_ids.items():
        summed_reps = feat_tab_reps[reps].sum(axis=1) # sum replicates for current sample
        final_tab_cols[sam] = summed_reps
    
    feat_tab_derep = pd.DataFrame(final_tab_cols)
    final_tab_df = pd.concat([feat_tab_derep, feat_tab_non_reps], axis=1) # join dereped sampled with singletons

    return final_tab_df

def map_reps_to_eblank(extr_blank_map):
    extr_blank_dict = defaultdict(list)
    with open(extr_blank_map) as f:
        f.readline()
        for line in f.readlines():
            split_line = line.strip().split("\t")
            if len(split_line) < 2: # if theres no entry in the eblank col, call it out
                sam_rep = line.strip()
                eblank_rep = "NA"
            else:
                sam_rep, eblank_rep = split_line
            extr_blank_dict[eblank_rep].append(sam_rep)

    return extr_blank_dict

def map_sams_to_fblank(meta_filepath_map, field_blank_map):
    sam_to_readfile_dict = {}
    with open(meta_filepath_map, "r") as f:
        f.readline()
        for line in f.readlines():
            sam_id, file_prefix = line.strip().split("\t")
            sam_to_readfile_dict[sam_id] = file_prefix.split(";")

    fblank_to_readfile_dict = defaultdict(list)
    with open(field_blank_map, "r") as f:
        f.readline()
        for line in f.readlines():
            sam_id, fblank = line.strip().split("\t")[:2]
            for p in sam_to_readfile_dict[sam_id]:
                if p == "NA":
                    continue
                else:
                    fblank_to_readfile_dict[fblank].append(p)

    return fblank_to_readfile_dict

def filter_by_abundance(feat_tab, n):
    sample_totals = feat_tab.sum(axis=0)
    feat_tab_filt_abun = feat_tab.div(sample_totals, axis=1)
    taxa_greater_than_n_abun = (feat_tab_filt_abun > n).sum(axis=1) > 0
    feat_tab_drop_low_abun = feat_tab[taxa_greater_than_n_abun]

    return feat_tab_drop_low_abun

def decontam(logger, feat_tab, study_metadata, extr_blank_map, meta_filepath_map, field_blank_map, contaminants, abun_cutoff, outdir):
    ftab_df = pd.read_csv(feat_tab, sep="\t", index_col="Taxon")

    extr_blank_dict = map_reps_to_eblank(extr_blank_map)
    ftab_minus_eblank = pd.DataFrame()
    no_eblank_reps = []
    for k, v in extr_blank_dict.items():
        try:
            subtracted_eblank = ftab_df[v].sub(ftab_df[k], axis=0)
        except KeyError:
            no_eblank_reps += v
            subtracted_eblank = ftab_df[v]
        ftab_minus_eblank = pd.concat([ftab_minus_eblank, subtracted_eblank], axis=1)

    ftab_eb_clipped = ftab_minus_eblank.clip(lower=0) # set values < 0 to 0
    ftab_eb_clipped.to_csv(outdir / "ftab_minus_eblank.tsv", sep="\t")

    ftab_derep_df = sum_replicates(ftab_df)
    ftab_derep_df.to_csv(outdir / "full_derep_ftab.tsv", sep="\t")

    field_blank_dict = map_sams_to_fblank(meta_filepath_map, field_blank_map)
    just_fblank_sams = pd.DataFrame()
    for fb in list(field_blank_dict.keys()):
        if fb in ftab_derep_df.columns:
            just_fblank_sams = pd.concat([just_fblank_sams, ftab_derep_df[fb]], axis=1)
    just_fblank_sams.index.rename("Taxon", inplace=True)
    just_fblank_sams.to_csv(outdir / "ftab_just_fblanks.tsv", sep="\t")

    ftab_minus_fblank = pd.DataFrame()
    no_fblank_sams = []
    for k, v in field_blank_dict.items():
        try:
            subtracted_fblank = ftab_derep_df[v].sub(ftab_derep_df[k], axis=0)
        except KeyError:
            no_fblank_sams += v
            subtracted_fblank = ftab_derep_df[v]
        ftab_minus_fblank = pd.concat([ftab_minus_fblank, subtracted_fblank], axis=1)

    ftab_fb_clipped = ftab_minus_fblank.clip(lower=0) # set values < 0 to 0
    ftab_fb_clipped.to_csv(outdir / "ftab_minus_fblank.tsv", sep="\t")

    ftab_drop_contams = ftab_fb_clipped.drop(contaminants, axis=0) # drop common contaminants

    zero_counts = ftab_drop_contams[ftab_drop_contams.sum(axis=1) < 1].index # get taxa that now have zero counts across all samples
    ftab_drop_zeros = ftab_drop_contams.drop(zero_counts, axis=0) # drop these taxa
    ftab_drop_zeros.to_csv(outdir / "ftab_no_contams.tsv", sep="\t")

    for n in abun_cutoff:
        feat_tab_drop_low_abun = filter_by_abundance(ftab_drop_contams, n)
        feat_tab_drop_low_abun.to_csv(outdir / f"feat_tab_{n}_abundance.tsv", sep="\t")

logger = "logger"
feat_tab = Path("input_decontam/DB_feat_tab.tsv").resolve()
study_metadata = Path("input_decontam/DB_metadata.tsv").resolve()
extr_blank_map = Path("input_decontam/DB_eblank_map.tsv").resolve()
meta_filepath_map = Path("input_decontam/DB_best_matches.tsv").resolve()
field_blank_map = Path("input_decontam/DB_field_blanks.tsv").resolve()
contaminants = [ # constricted based on dambaseline output
        "Homo sapiens", # humans
        "Bos taurus", # domestic cow
        "Gallus gallus", # domestic chicken
        "Squatina", # Angelshark
        "Felis catus", # domestic cat
] # left in Canis, Canidae (could be wolves, dogs, idk)
abun_cutoff = [0.001, 0.005, 0.01, 0.05, 0.1]

outdir = Path("test_decontam_out")
outdir.mkdir(exist_ok=True)
decontam(logger, feat_tab, study_metadata, extr_blank_map, meta_filepath_map, field_blank_map, contaminants, abun_cutoff, outdir)