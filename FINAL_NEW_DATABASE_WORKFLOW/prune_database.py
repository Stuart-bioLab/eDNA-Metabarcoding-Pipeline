import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from collections import defaultdict
import re

def read_file_data(in_seq_file, in_tax_file):
    """Read in data from database files and concatenate taxonomy with sequence."""
    tax_df = pd.read_csv(in_tax_file, sep="\t", index_col="Feature ID")
    tax_df["Taxon"] = tax_df["Taxon"].apply(lambda x: re.sub(r"k__.*?;", "", x)) # kingdoms are NA sometimes
    dict_taxa = tax_df["Taxon"].to_dict()

    species_names = tax_df["Taxon"].apply(
        lambda x: x.split(";g__")[-1].replace(";s__", " ").replace("_", " ")
        ).sort_values().unique().tolist()
    with open("gbif_species_hits.txt", "w") as f:
        for sp in species_names:
            f.write(f"{sp}\n")

    records = SeqIO.parse(in_seq_file, format="fasta")
    for r in records:
        dict_taxa[r.name] += "__seq__" + str(r.seq)

    return dict_taxa

def add_species(dict_taxa, missing_sp, inhs_seq_file, inhs_tax_file):
    """Add species from the full INHS database if they're in the Fox but were missed by GBIF."""
    feats_to_add = []

    with open(inhs_tax_file) as f:
        for line in f.readlines():
            if any(sp in line for sp in missing_sp):
                feat, tax = line.split("\t")
                feats_to_add.append(feat)
                dict_taxa[feat] = re.sub(r"k__.*?;", "", tax.strip())
    
    with open(inhs_seq_file) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            feat = lines[i].replace(">", "").strip()
            if feat in missing_sp:
                dict_taxa[feat] += "__seq__" + lines[i+1].strip()
    
    return dict_taxa

def prune_files(dict_taxa):
    """For entries that have identical taxonomy and sequence, keep only one."""
    pruned_dict = {}

    for k, v in dict_taxa.items():
        if v not in pruned_dict.values(): # only add entry if theres no other entry with identical taxonomy and sequence
            pruned_dict[k] = v

    pruned_df = pd.DataFrame(pruned_dict.items(), columns=["Feature ID", "Taxon"])
    pruned_df["Taxon"] = pruned_df["Taxon"].apply(lambda x: "k__Metazoa;" + x.split("__seq__")[0]) # put kingdom assignment back in
    pruned_df.to_csv("gbif_pruned_tax.tsv", sep="\t", index=False)

    with open("gbif_pruned_seq.fasta", "w") as f:
        for k, v in pruned_dict.items():
            seq = v.split("__seq__")[-1]
            f.write(f">{k}\n{seq}\n")
    
    return

def main():
    in_seq_file = "gbif_seq.fasta"
    in_tax_file = "gbif_tax.tsv"

    inhs_seq_file = "../updated_12S_database_files/dna-sequences.fasta"
    inhs_tax_file = "../updated_12S_database_files/inhs_plus_picq_stringent.tsv"

    dict_taxa = read_file_data(in_seq_file, in_tax_file)

    missing_sp = ["molitrix"]
    dict_taxa = add_species(dict_taxa, missing_sp, inhs_seq_file, inhs_tax_file)
    
    prune_files(dict_taxa)

def read_inhs_files(seq_file, tax_file):
    """Read in data from database files and concatenate taxonomy with sequence."""
    tax_df = pd.read_csv(tax_file, sep="\t", index_col="Feature ID")
    tax_df["Taxon"] = tax_df["Taxon"].apply(lambda x: re.sub(r"k__.*?;", "", x)) # kingdoms are NA sometimes
    dict_taxa = tax_df["Taxon"].to_dict()

    species_names = tax_df["Taxon"].apply(
        lambda x: x.split(";g__")[-1].replace(";s__", " ").replace("_", " ")
        ).sort_values().unique().tolist()
    with open("gbif_species_hits.txt", "w") as f:
        for sp in species_names:
            f.write(f"{sp}\n")

    records = SeqIO.parse(seq_file, format="fasta")
    for r in records:
        dict_taxa[r.name] += "__seq__" + str(r.seq)

    return dict_taxa

def main():
    inhs_seq = "../updated_12S_database_files/dna-sequences.fasta"
    inhs_tax = "../updated_12S_database_files/taxonomy.tsv"
    species_list = "FINAL_gbif_add_remove_metzke_refs.txt"

    read_inhs_files(inhs_seq, inhs_tax)

main()