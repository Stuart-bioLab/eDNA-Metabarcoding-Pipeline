# Update picq database with space-separated list of genbank ids

import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from Bio import Entrez, SeqIO
import argparse
from collections import defaultdict

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tax_file")
    parser.add_argument("-s", "--seq_file")
    parser.add_argument("-i", "--id_list", nargs="+")
    args = parser.parse_args()
    return args

def fetch_seqs(ids):
    new_entries = defaultdict(lambda: {"seq": None, "tax_id": None, "tax_str": None})

    handle = Entrez.efetch(
        db="nucleotide",
        id=ids,
        rettype="gb",
        retmode="text"
    )

    records = list(SeqIO.parse(handle, "genbank"))
    for r in records:
        new_entries[r.id]["seq"] = r.seq
        for feat in r.features: # record features store genetic and taxonomic info
            if feat.type == "source": # looking for taxonomy information
                new_entries[r.id]["tax_id"] = feat.qualifiers.get("db_xref")[0][6:] # just get the taxonomy id number
    
    return new_entries

def fetch_tax(new_entries):
    tax_id_list = [new_entries[x]["tax_id"] for x in new_entries] # extract all taxonomic ids
    handle = Entrez.efetch(
        db="taxonomy",
        id=tax_id_list,
        retmode="xml"
    )

    record = Entrez.read(handle)

    tax_id_map = {}
    for entry in record:
        tax_id = entry.get("TaxId") # get the id to map back to the hits dict
        tax_dict = { # initialize dict to store taxonomy info
            "kingdom": "k__NA",
            "phylum": "p__NA",
            "class": "c__NA",
            "order": "o__NA",
            "family": "f__NA",
            "genus": "g__NA",
            "species": "s__NA"
        }
        for item in entry.get("LineageEx", []): # extract tax info and format for qiime
            rank = item["Rank"]
            name = item["ScientificName"]
            if rank in tax_dict.keys():
                tax_dict[rank] = f"{rank[0]}__{name}"
        species = entry.get("ScientificName", []).split(" ")[1]
        tax_dict["species"] = "s__" + species
        tax_id_map[tax_id] = tax_dict

    return tax_id_map

def export_database(seq_file, tax_file, outdir):
    """Decompress qiime archive files so that the database files can be written to."""
    subprocess.run([ # decompress sequences
        "qiime", "tools", "export",
        "--input-path", seq_file,
        "--output-path", outdir
    ])

    subprocess.run([ # decompress taxonomy
        "qiime", "tools", "export",
        "--input-path", tax_file,
        "--output-path", outdir
    ])

def append_entries(db_dict, seq_file, tax_file):
    """Add sequence and taxonomy for input records to respective file."""
    with open(seq_file, "a") as f:
        for record_id in db_dict:
            seq = db_dict[record_id]["seq"]
            f.write(f">{record_id}\n{seq}\n")
    
    with open(tax_file, "a") as f:
        for record_id in db_dict:
            tax_str = db_dict[record_id]["tax_str"]
            joined_tax_str = ";".join(list(tax_str.values()))
            f.write(f"{record_id}\t{joined_tax_str}\n")

def import_database(seq_file, tax_file, outdir):
    """Compress sequence and taxonomy files back to qiime archive format."""
    subprocess.run([ # compress sequences
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", seq_file,
        "--output-path", outdir / "updated_seq"
    ])

    subprocess.run([ # compress taxonomy
        "qiime", "tools", "import",
        "--type", "FeatureData[Taxonomy]",
        "--input-path", tax_file,
        "--output-path", outdir / "updated_tax"
    ])

def main():
    args = get_args()
    old_tax = args.tax_file
    old_seq = args.seq_file
    ids = args.id_list
    outdir = Path("updated_database").resolve()
    outdir.mkdir(exist_ok=True)

    Entrez.email = "bmoginot5@gmail.com"

    new_entries = fetch_seqs(ids)
    tax_id_map = fetch_tax(new_entries)
    
    for record_id in new_entries.keys():
        new_entries[record_id]["tax_str"] = tax_id_map[new_entries[record_id]["tax_id"]]

    export_database(old_seq, old_tax, outdir)
    unzipped_seq = outdir / "dna-sequences.fasta"
    unzipped_tax = outdir / "taxonomy.tsv"
    append_entries(new_entries, unzipped_seq, unzipped_tax)
    import_database(unzipped_seq, unzipped_tax, outdir)

if __name__ == "__main__":
    main()
