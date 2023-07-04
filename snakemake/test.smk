#config
configfile: "./snakemake/config.yaml"


FILES = json.load(open(config['SAMPLES_JSON']))

print(FILES)

import csv
import os

SAMPLES = sorted(FILES.keys())
print(SAMPLES)
# List all samples by sample_name, sample_type, and assay
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type+ "_" + assay)

print(MARK_SAMPLES)

CUT_TAG = config["c_t"]
CHIP = config["chip"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP not in sample]
import json
print(CUT_TAGS)

fastq_paths = FILES[sample][sample_type][assay]

print(fastq_paths)
# Rule without wildcards at the top of the workflow
rule all:
    input:
        expand("cat/{sample}_R1.fastq", sample=CUT_TAGS),
        expand("cat/{sample}_R2.fastq", sample=CUT_TAGS)

rule merge_fastq:
    input:
        R1 = [path for path in fastq_paths if "_R1_" in path],
        R2 = [path for path in fastq_paths if "_R1_" in path]
    output:
        "cat/{sample}_R1.fastq",
        "cat/{sample}_R2.fastq"
    threads: 20
    shell:
        """
        gunzip -c {input.R1} | gzip > {output[0]} 
        gunzip -c {input.R2} | gzip > {output[1]} 
        """