# Remove or append species to picq database from list

import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import argparse
from collections import defaultdict
import edlib
import pandas as pd
from Bio.SeqIO.FastaIO import FastaWriter

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tax_file")
    parser.add_argument("-s", "--seq_file")
    parser.add_argument("-n", "--names", nargs="+")
    parser.add_argument("-e", "--email")
    parser.add_argument("-a", "--append", action="store_true")
    parser.add_argument("-r", "--remove", action="store_true")
    parser.add_argument("-f", "--file")
    args = parser.parse_args()
    return args

def parse_names(names_file):
    """Get list of names from input file."""
    with open(names_file, "r") as f:
        names_list = [x.strip() for x in f.readlines()]

    return names_list

def search_name(name):
    """Seach for mitochondrial sequences associated with target taxon and get first 20 results."""
    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{name}[organism] AND rrn12[gene]",
        retmax=20
    )

    results = Entrez.read(handle)
    return results["IdList"]

def fetch_seq(name, ids, hits):
    """Iterate over IDs until a 12S sequence is found. Store record ID, taxonomy ID, and sequence"""
    handle = Entrez.efetch(
        db="nucleotide",
        id=ids,
        rettype="gb",
        retmode="text"
    )

    records = list(SeqIO.parse(handle, "genbank"))
    for r in records:
        for f in r.features:
            if f.type == "gene": # only looking for mitchondrial rRNA sequences
                gene_name = f.qualifiers.get("gene", []) # get the name of the gene
                if "rrn12" in gene_name: # only look at record if its a 12S sequence
                    hits[name]["seq"] = f.extract(r.seq) # store 12S sequence
                    hits[name]["record_id"] = r.id
                    hits[name]["label"] = f"{name} {r.id}".replace(" ", "_")
            if f.type == "source": # looking for taxonomy information
                hits[name]["tax_id"] = f.qualifiers.get("db_xref")[0][6:] # just get the taxonomy id number
        if hits[name]["seq"]: # if 12S sequence found, exit loop
            print(f"Found 12S sequence for {name}.")
            break
    else:
        print(f"Could not find sequence for {name}.")
        del hits[name]

    return hits

def find_v5(f_primer, r_primer, target_seq):
    """Align primers to 12S sequence and slice out substring between them."""
    result = edlib.align(f_primer, target_seq, mode="HW", task="path")
    start = result["locations"][0][1] + 1

    result = edlib.align(r_primer, target_seq, mode="HW", task="path")
    end = result["locations"][0][0]

    return target_seq[start:end]

def fetch_tax(hits):
    """Get taxonomy information for extracted record."""
    name_map = {hits[x]["tax_id"]: x for x in hits} # store name with tax id so its easy to map back to hits dict
    tax_id_list = [hits[x]["tax_id"] for x in hits] # extract all taxonomic ids
    handle = Entrez.efetch(
        db="taxonomy",
        id=tax_id_list,
        retmode="xml"
    )

    records = Entrez.read(handle)

    for r in records:
        tax_id = r.get("TaxId") # get the id to map back to the hits dict
        tax_dict = { # initialize dict to store taxonomy info
            "kingdom": "k__NA",
            "phylum": "p__NA",
            "class": "c__NA",
            "order": "o__NA",
            "family": "f__NA",
            "genus": "g__NA",
            "species": "s__NA"
        }
        for item in r.get("LineageEx", []): # extract tax info and format for qiime
            rank = item["Rank"]
            name = item["ScientificName"]
            if rank in tax_dict.keys():
                tax_dict[rank] = f"{rank[0]}__{name}"
        species = r.get("ScientificName", []).split(" ")[1]
        tax_dict["species"] = "s__" + species

        name = name_map[tax_id]
        hits[name]["tax_str"] = ";".join(tax_dict.values())

        print(f"Found taxonomy for {name}.")

    return hits

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

def append_entries(hits, seq_file, tax_file, f_primer, r_primer):
    """Add sequence and taxonomy for input records to respective file."""
    with open(seq_file, "a") as f:
        for name in hits:
            label = hits[name]["label"]
            seq = hits[name]["seq"]
            f.write(f">{label}\n{seq}\n")
    
    with open(tax_file, "a") as f:
        for name in hits:
            label = hits[name]["label"]
            tax_str = hits[name]["tax_str"]
            f.write(f"{label}\t{tax_str}\n")

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

def append(email, names_to_search, unzipped_seq, unzipped_tax):
    Entrez.email = email
    if not Entrez.email:
        print("Please supply Entrez with an email using --email argument.")
        sys.exit(1)

    hits = defaultdict(lambda: {"record_id": None, "seq": None, "tax_id": None, "tax_str": None, "label": None})

    for name in names_to_search:
        print(f"Searching for {name}...")
        target_ids = search_name(name)
        hits = fetch_seq(name, target_ids, hits)
        if name in hits.keys():
            print(f"Finding V5 region for {name}...")
            hits[name]["seq"] = find_v5(f_primer, r_primer, hits[name]["seq"])
    
    print("Finding taxonomy for all entries...")
    hits = fetch_tax(hits)

    print("Appending sequences and taxonomy...")
    append_entries(hits, unzipped_seq, unzipped_tax, f_primer, r_primer)

def remove(names_to_search, unzipped_seq, unzipped_tax, outdir):
    tax_df = pd.read_csv(unzipped_tax, sep="\t")
    tax_df["sp_name"] = tax_df["Taxon"].apply(lambda x: x.split(";g__")[-1].replace(";s__", " "))
    target_rows = tax_df["sp_name"].isin(names_to_search)

    pruned_tax_file = tax_df[~target_rows]
    pruned_tax_file[["Feature ID", "Taxon"]].to_csv( outdir / "pruned_tax_file.tsv", sep="\t", index=False)

    dropped_rows = tax_df[target_rows]
    target_labels = dropped_rows["Feature ID"].to_list()

    parsed_fasta = SeqIO.parse(unzipped_seq, format="fasta")
    pruned_records = [record for record in parsed_fasta if record.name not in target_labels]

    with open(outdir / "pruned_seq_file.fasta", "w") as f:
        writer = FastaWriter(f, wrap=None)
        writer.write_file(pruned_records)

    return

def main():
    args = get_args()

    old_tax = args.tax_file
    if not old_tax:
        old_tax = Path("picq_2024_12S_database/db_tax.qza")

    old_seq = args.seq_file
    if not old_seq:
        old_seq = Path("picq_2024_12S_database/db_seq.qza")

    names_to_search = args.names
    if not names_to_search:
        if not args.file:
            print("Please give taxa to search for with either of the --names or --file arguments.")
            sys.exit(1)
        infile = Path(args.file)
        if not infile.exists():
            print("Could not find file.")
            sys.exit(1)
        names_to_search = parse_names(infile)

    f_primer = Seq("ACTGGGATTAGATACCCC")
    r_primer = Seq("TAGAACAGGCTCCTCTAG").reverse_complement()

    outdir = Path("updated_database").resolve()
    outdir.mkdir(exist_ok=True)

    # print("Decompressing input files...")
    # export_database(old_seq, old_tax, outdir)
    unzipped_seq = outdir / "dna-sequences.fasta"
    unzipped_tax = outdir / "taxonomy.tsv"

    if args.append:
        append(args.email, names_to_search, unzipped_seq, unzipped_tax)
    
    if args.remove:
        remove(names_to_search, unzipped_seq, unzipped_tax, outdir)

    # print("Re-importing updated database...")
    # import_database(unzipped_seq, unzipped_tax, outdir)

    print("All done.")

if __name__ == "__main__":
    main()
