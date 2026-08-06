# Fox-River-eDNA-metabarcoding
An eDNA metabarcoding project focused on identifying and analyzing biodiversity through environmental DNA sequencing. This workflow uses molecular markers and bioinformatics approaches to detect species from environmental samples without direct organism collection.
The pipeline architecture is based on those built by B. Moginot (https://github.com/bmoginot/Fox-River-eDNA-Pipeline) and R. Patel (https://github.com/richapatel138/WildMileeDNAPipeline).

# ETL framework
Data are in disparate folders, not organized by study. The pipeline needs a manifest which lists the paths to read files in a text file. The metadata maps each sample to a study(s), so it can be used to figure out where the read files for the samples are located.  
## Pre-processing
The script `make_manifest_metadata.py` takes the metadata, study type, and read data directory as input. It subsets the metadata for samples from the target study. It then searches the input directory and pulls out paths to files that are labelled with the sample ids from the target study. The script then creates a manifest for the pipeline. It then uses the manifest ids to index the metadata, pairing replicates with samples and generating new metadata for just the target samples. The new metadata and manifest files are required as input for the pipeline.  
`make_manifest_metadata.py` handles duplicate values in the metadata as well as labelled duplicate read files. *This data handling is specific to this project and this metadata.*  

This may be implemented into the pipeline so that only one script needs to be run, but we'll cross that bridge when we get there.
