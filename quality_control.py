import subprocess
import sys

manifests = sys.argv[1:]
# manifests = ["090326_runs/DamBaseline_full_run_with_MHD-DMN/final_DamBaseline_manifest.tsv", "090426_runs/JuneJulyTemporal_full/final_JuneJulyTemporal_manifest.tsv", "090426_runs/EbonyTemporal_full/final_EbonyTemporal_manifest.tsv", "090426_runs/Filter_full_run/final_Filter_5.0v0.45_manifest.tsv"]

final_manif = "all_studies_manifest.tsv"

with open(final_manif, "w") as f:
    f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")

for manif in manifests:
    with open(final_manif, "r") as f:
        written_reads = f.readlines()

    with open(manif, "r") as m:
        manif_lines = m.readlines()

    with open(final_manif, "a") as f:
        for line in manif_lines:
            if line not in written_reads:
                f.write(line)

reads = []
with open(final_manif, "r") as f:
    f.readline()
    for line in f.readlines():
        sam_id, fpath, rpath = line.strip().split("\t")
        reads.append(fpath)
        reads.append(rpath)

for r in reads:
    subprocess.run([
        "fastqc",
        "--outdir", "fastqc_reports",
        "--threads", "20",
        r
    ])

subprocess.run([
    "multiqc", "fastqc_reports"
])