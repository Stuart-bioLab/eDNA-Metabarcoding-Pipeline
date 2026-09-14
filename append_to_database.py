# Update picq database with newly-generated mitochondrial sequences

import subprocess
import sys
from collections import defaultdict
from pathlib import Path

def export_database(seq_file, tax_file, db_dir):
    subprocess.run([
        "qiime", "tools", "export",
        "--input-path", seq_file,
        "--output-path", db_dir
    ])

    subprocess.run([
        "qiime", "tools", "export",
        "--input-path", tax_file,
        "--output-path", db_dir
    ])

def read_entries(seq_file, tax_file):
    db_dict = defaultdict(lambda: {"seq": None, "tax": None})
    with open(seq_file, "r") as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            taxid = lines[i][1:]
            sequence = lines[i+1]
            seq_dict[taxid]["seq"] = sequence
    
    with open(tax_file, "r") as f:
        f.readline()
        for line in f.readlines():
            taxid, taxonomy = line.strip().split("\t")
            db_dict[taxid]["tax"] = taxonomy
    
    return db_dict

def append_entries(db_dict, seq_file, tax_file):
    with open(seq_file, "a") as f:
        for taxid in db_dict:
            f.write(f">{taxid}\n{db_dict[taxid]["seq"]}\n")
    
    with open(tax_file, "a") as f:
        for taxid in db_dict:
            f.write(f"{taxid}\t{db_dict[taxid]["tax"]}\n")

def import_database(seq_file, tax_file, outdir):
    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", seq_file,
        "--output-path", outdir / "updated_picq_db_seq"
    ])

    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Taxonomy]",
        "--input-path", tax_file,
        "--output-path", outdir / "updated_picq_db_tax"
    ])

def main():
    picq_db_dir = Path("pdb_test").resolve()
    picq_db_seq = picq_db_dir / "db_seq.qza"
    picq_db_tax = picq_db_dir / "db_tax.qza"

    export_database(picq_db_seq, picq_db_tax, picq_db_dir)
    db_dict = read_entries(sys.argv[1], sys.argv[2])
    unzipped_seq = picq_db_dir / "dna-sequences.fasta"
    unzipped_tax = picq_db_dir / "taxonomy.tsv"
    append_entries(db_dict, unzipped_seq, unzipped_tax)
    import_database(unzipped_seq, unzipped_tax, picq_db_dir)

if __name__ == "__main__":
    main()
