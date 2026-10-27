# Update picq database with space-separated list of species names

import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import argparse
from collections import defaultdict

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tax_file")
    parser.add_argument("-s", "--seq_file")
    parser.add_argument("-n", "--names", nargs="+")
    args = parser.parse_args()
    return args

def search_name(name):
    """Seach for mitochondrial sequences associated with target taxon and get first 20 results."""
    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{name}[organism] AND rrn12[gene]",
        retmax=20
    )

    results = Entrez.read(handle)
    return results["IdList"]

def fetch_seqs(name, ids, hits):
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
            print(hits[name])
            label = hits[name]["label"]
            seq = hits[name]["seq"]
            v5_region = seq.split(f_primer)[-1].split(r_primer)[0] # this is really loosey-goosey
            f.write(f">{label}\n{v5_region}\n")
    
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

def main():
    args = get_args()

    old_tax = args.tax_file
    if not old_tax:
        old_tax = "picq_2024_12S_database/db_tax.qza"

    old_seq = args.seq_file
    if not old_seq:
        old_seq = "picq_2024_12S_database/db_seq.qza"

    names_to_search = args.names
    if not names_to_search:
        print("Please give taxa to search for with --names argument.")
        sys.exit(1)

    f_primer = Seq("ACTGGGATTAGATACCCC")
    r_primer = Seq("TAGAACAGGCTCCTCTAG").reverse_complement()[:-3]

    outdir = Path("updated_database").resolve()
    outdir.mkdir(exist_ok=True)

    Entrez.email = "bmoginot5@gmail.com"

    hits = defaultdict(lambda: {"record_id": None, "seq": None, "tax_id": None, "tax_str": None, "label": None})

    for name in names_to_search:
        print(f"Searching for {name}...")
        target_ids = search_name(name)
        hits = fetch_seqs(name, target_ids, hits)
    
    print("Finding taxonomy for all entries...")
    hits = fetch_tax(hits)

    # print("Decompressing input files...")
    # export_database(old_seq, old_tax, outdir)
    unzipped_seq = outdir / "dna-sequences.fasta"
    unzipped_tax = outdir / "taxonomy.tsv"

    print("Appending sequences and taxonomy...")
    append_entries(hits, unzipped_seq, unzipped_tax, f_primer, r_primer)

    # print("Re-importing updated database...")
    # import_database(unzipped_seq, unzipped_tax, outdir)

    print("All done.")

if __name__ == "__main__":
    main()
