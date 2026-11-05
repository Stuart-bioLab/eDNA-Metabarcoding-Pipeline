# CRABS commands used to retrieve missing Illinois seqs and add them to the database

## Set up crabs
### Install executable
cloned repo from github into Software/ folder  
made a symlink to the crabs executable via  
```
sudo ln -s /home/fishl/Software/reference_database_creator/crabs /usr/local/bin
```

### Install dependencies
created new env and installed dependencies via  
```
conda create -n crabs -c conda-forge -c bioconda \
python=3.11 requests rich rich-click matplotlib numpy \
blast cutadapt vsearch clustalw fasttree
```

### Store files in this dir
```
mkdir new_crabs_db
cd new_crabs_db
```

## Download databases
### Download NCBI taxonomy
```
crabs --download-taxonomy --output .
```

### Download NCBI sequences
grab vertebrate sequences containing the 12S gene from ncbi.
```
crabs --download-ncbi \
  --database nucleotide \
  --query '(Vertebrata[Organism] AND mitochondrion[filter] AND ("12S ribosomal RNA"[Title] OR 12S[Gene] OR rrnS[All Fields]) AND 100:25000[SLEN])' \
  --email bmoginot5@gmail.com \
  --output ncbi_12S.fasta
```

import to crabs format
```
crabs --import \
  --import-format ncbi \
  --input ncbi_12S.fasta \
  --names names.dmp \
  --nodes nodes.dmp \
  --acc2tax nucl_gb.accession2taxid \
  --output ncbi_12S_crabs.txt
```

### Download MIDORI sequences
crabs can't access midori's website for some reason. not even when i try the version opale downloaded. instead, were going to download the raw fasta from zenodo
```
wget "https://zenodo.org/records/22685701/files/MIDORI2_GB273_RAW_sp.zip?download=1"
```

unzip the 12s fasta file
```
zip -dc RAW_sp/uniq/MIDORI2_UNIQ_SP_NUC_GB273_srRNA_RAW.fasta.gz > midori_12S.fasta
```

import the fasta file into crabs format
```
crabs --import \
  --import-format MIDORI \
  --input midori_12S.fasta \
  --names names.dmp \
  --nodes nodes.dmp \
  --acc2tax nucl_gb.accession2taxid \
  --output midori_12S_crabs.txt \
  --ranks 'superkingdom;phylum;class;order;family;genus;species'
```

### Download MitoFish sequences
download with the native crabs function
```
crabs --download-mitofish --output mitofish.fasta
```

import into crabs format
```
crabs --import \
  --import-format mitofish \
  --input mitofish.fasta \
  --names names.dmp \
  --nodes nodes.dmp \
  --acc2tax nucl_gb.accession2taxid \
  --output mitofish_crabs.txt
```

## Merge and Filter databases
### Merge downloaded seqs into one file
```
crabs --merge \
  --input 'ncbi_12S_crabs.txt;midori_12S_crabs.txt;mitofish_crabs.txt' \
  --uniq \
  --output merged_crabs.txt
```

### In-silico PCR
extract V5 region using Riaz primers
```
crabs --in-silico-pcr \
  --input merged_crabs.txt \
  --output insilico_pcr.txt \
  --forward ACTGGGATTAGATACCCC \
  --reverse TAGAACAGGCTCCTCTAG
```

recover sequences using pairwise global alignment
```
crabs --pairwise-global-alignment \
  --input merged_crabs.txt \
  --amplicons insilico_pcr.txt \
  --output aligned_crabs.txt \
  --forward ACTGGGATTAGATACCCC \
  --reverse TAGAACAGGCTCCTCTAG \
  --percent-identity 0.95 \
  --coverage 95
```

### Derepelicate
remove duplicate sequences
```
crabs --dereplicate \
  --input aligned_crabs.txt \
  --output dereplicated_crabs.txt \
  --dereplication-method unique_species
```

### Subset
check which species are actually in the crabs database
```
zero_count=0

while read species; do
    count=$(grep -Fc "$species" aligned_crabs.txt)
    if [ "$count" -gt 0 ]; then
        echo "$species"
        ((zero_count++))
    fi
done < ../database_files_103026/idnr_set.txt

echo "Species found: $zero_count"
```

subset based on species list from IDNR
```
crabs --subset \
  --input aligned_crabs.txt \
  --output target_species_crabs.txt \
  --include ../database_files_103026/idnr_set.txt
```

### Filter out ambiguous seqs
```
crabs --filter \
  --input target_species_crabs.txt \
  --output filtered_crabs.txt \
  --maximum-n 10
```

## Convert to QIIME format
### Export
export sequence file
```
crabs --export \
  --input filtered_crabs.txt \
  --output qiime_reference.fasta \
  --export-format qiime-fasta
```

export taxonomy file
```
crabs --export \
  --input filtered_crabs.txt \
  --output qiime_taxonomy.txt \
  --export-format qiime-text
```

format taxonomy file
```
sed "s/;s__.*_/;s__/g" qiime_taxonomy.txt > edited_tax.tsv
```

### Append to existing database
combine sequences with those in existing seq file
```
cat ../database_files_103026/dna-sequences.fasta qiime_reference.fasta > final_seq.fasta
```

combine taxonomy with existing tax file
```
cat ../database_files_103026/taxonomy.tsv edited_tax.tsv > final_tax.tsv
```

### Import as QIIME archive
import seqs
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path final_seq.fasta \
  --output-path final_seq.qza
```

import tax
```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path final_tax.tsv \
  --output-path final_tax.qza
```