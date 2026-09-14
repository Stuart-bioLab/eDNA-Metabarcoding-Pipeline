# Fox-River-eDNA-metabarcoding
An eDNA metabarcoding project focused on identifying and analyzing biodiversity through environmental DNA sequencing. This workflow uses molecular markers and bioinformatics approaches to detect species from environmental samples without direct organism collection.
This pipeline builds on the work of [B. Moginot](https://github.com/bmoginot/Fox-River-eDNA-Pipeline) and [R. Patel](https://github.com/richapatel138/WildMileeDNAPipeline) and is ultimately based on the analysis done by [Picq et al. 2024](https://doi.org/10.1080/02705060.2024.2382454).

# Running this Pipeline (in General)
Detailed steps for running this pipeline can be found in the wiki, but here is a general overview.
## Locate Input Files
To run this pipeline you will need the following. More information on input can be found on the [Data]() page.
1. A directory containing sequence read files for the study of interest. These are likely contained in an external drive which should be mounted to the PC. See [Data]() for information on storage and access of this data.
2. A Sample Metadata file which lists information about collected samples and their related sequence reads. This file is likely found in your current working directory.
3. An Extraction Blank Table file that pairs each sample with a corresponding extraction blank. This file is likely found in your current working directory.
## Generate Preliminary Files
At the time of writing this README, not all samples for this project have been sequenced. As sequencing progresses, certain input files will need to be updated to facilitate the processing of new reads. Running the `locate_samples.py` script creates up to date input files that include all currently sequenced samples. Namely, this script generates a Manifest and a Blank Metadata file, which are required as input for the pipeline. Refer to [Pre-processing](https://github.com/Stuart-bioLab/eDNA-Metabarcoding-Pipeline/wiki/Pre%E2%80%90processing) for more information.
## Optional Steps
You may also want to update databases with new sequences or perform quality control before running the pipeline. These steps are not required, but they are detailed in [Pre-processing](https://github.com/Stuart-bioLab/eDNA-Metabarcoding-Pipeline/wiki/Pre%E2%80%90processing).
## Running Analysis
Once the **Study Manifest** and **Blank Metadata** files have been generated, they can be passed to `metabarcoding_pipeline.py` with the `--manifest` and `--blank_metadata` arguments, respectively. The remaining arguments will be handled by the **config.ini** file. Refer to [Data]() for more information on the config file.
## View Output Files
Final feature tables will be found in `results/final_output/`. Otherwise, there are several intermediate output files in the subdirectories of results/ which may need to be viewed to diagnose pipeline issues. Refer to [Data]() for more information on output files.
