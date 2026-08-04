# Fox-River-eDNA-metabarcoding
An eDNA metabarcoding project focused on identifying and analyzing biodiversity through environmental DNA sequencing. This workflow uses molecular markers and bioinformatics approaches to detect species from environmental samples without direct organism collection.

# ETL framework
Data are in disparate folders, not organized by study. The pipeline needs a manifest which lists the paths to read files in a text file. The metadata maps each sample to a study(s), so it can be used to figure out where the read files for the samples are located.
The script `find_data.py` takes the metadata, study type, and read data directory as input. It subsets the metadata for samples from the target study. It then searches the input directory and pulls out paths to files that are labelled with the sample ids from the target study.
The script then creates a manifest for the pipeline.
*This might be implemented in the pipeline with an option to give just a data directory*.
