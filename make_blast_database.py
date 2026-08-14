# Retrieve all 12S vertebrate sequences from NCBI, create database files for the qiime pipeline

from Bio import Entrez, SeqIO
import Bio

def generate_fasta(batch):
    """Fetch all records that matched query and filter for 12S sequences"""
    handle = Entrez.efetch(
        db="nucleotide",
        id=batch,
        rettype="gb",
        retmode="text"
    )

    records = list(SeqIO.parse(handle, "genbank"))

    hits = {} # dict of dicts to attribute sequence and taxonomy info to record id
    maxhits = len(batch)
    i = 0
    j = 0
    k = 0
    for record in records: # iterate over fetched entries
        if not record.seq.defined: # only look at record if it has a defined sequence
            k += 1
        else:
            if record.id not in hits.keys():
                hits[record.id] = {} # create entry for record
            for feat in record.features: # record features store genetic and taxonomic info
                if feat.type == "source": # looking for taxonomy information
                    tax_id = feat.qualifiers.get("db_xref")[0][6:] # just get the taxonomy id number
                    hits[record.id]["tax_id"] = tax_id # get the tax id and store it with the record id
                    j += 1
                if feat.type == "rRNA": # only looking for mitchondrial rRNA sequences
                    product = feat.qualifiers.get("product", []) # product stores type of rRNA gene
                    matches = ["12S ribosomal RNA", "s-rRNA", "rrnS"] # only looking for 12S sequences
                    if any(x in product for x in matches): # only look at record if its a 12S sequence
                        hits[record.id]["seq"] = feat.extract(record.seq) # store the 12S sequence for this record
                        i += 1
    print("number of 12S seqs:", i, "/", maxhits)
    print("number of tax ref ids:", j, "/", maxhits)
    print("number of undefined seqs:", k, "/", maxhits)
    print("number of entries created:", len(hits))

    return hits

def generate_taxonomy(hits):
    """Fetch taxonomy information for retrieved records and store in qiime format"""
    tax_id_list = [hits[x]["tax_id"] for x in hits] # extract all taxonomic ids
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
        tax_dict["kingdom"] = "k__NA"
        genus, species = entry.get("ScientificName", []).split(" ")[:2]
        tax_dict["genus"] = "g__" + genus
        tax_dict["species"] = "s__" + genus + "_" + species
        tax_id_map[tax_id] = tax_dict
    
    return tax_id_map

def write_files(hits, tax_map):
    """Write out sequence and taxonomy information to fasta and qiime taxa format, respectively."""
    with open("blast_12S_db_seq.fasta", "a") as f:
        for record in hits:
            seq = hits[record]["seq"]
            f.write(f">{record}\n")
            f.write(f"{seq}\n")
    
    with open("blast_12S_db_tax.tsv", "a") as t:
        for record in hits:
            tax_id = hits[record]["tax_id"]
            tax_str = tax_map[tax_id]
            joined_tax_str = ";".join(list(tax_str.values()))
            t.write(f"{record}\t{joined_tax_str}\n")

def main():
    query = "refseq[filter] AND vertebrates[organism] AND mitochondrion[filter]"

    Entrez.email = "bmoginot@gmail.com"

    handle = Entrez.esearch(
        db="nucleotide",
        term=query,
        retmax=100000
    )

    results = Entrez.read(handle)

    print("total hits:", results["Count"])

    ids = results["IdList"]

    batch_size = 500
    j = 1
    for i in range(0, len(ids), batch_size): # batch queries to not be booted from NCBI
        print("batch number:", j)
        batch = ids[i:i + batch_size]
        j += 1
    
        hits = generate_fasta(batch)

        hits = { # failsafe to remove any entries that did not have a sequence
            record: data
            for record, data in hits.items()
            if "seq" in data
        }
        print(f"number of entries after cleaning:", len(hits))

        tax_map = generate_taxonomy(hits)

        write_files(hits, tax_map)
    print("FINISHED")

main()