########################################
# Load config and sample definitions
########################################
import json
import os

configfile: "/home/mattia/cfChromatin_in_MS/snakemake/config_ChIP.yaml"

# Load JSON with sample definitions
FILES = json.load(open(config['SAMPLES_JSON']))

# Sort samples
SAMPLES = sorted(FILES.keys())

# Build list of all sample+mark combinations
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type + "_" + assay)

# Filter by assay type
CUT_TAG   = config["c_t"]
CHIP      = config["chip"]

CUT_TAGS  = [s for s in MARK_SAMPLES if CUT_TAG in s]
CHIPS     = [s for s in MARK_SAMPLES if CHIP in s]
CHIPS_SE  = [s for s in MARK_SAMPLES if CHIP in s]  # adjust if needed

RUNID = config["RUN_ID"]

########################################
# Build target file lists
########################################

ALL_SAMPLES = CHIPS + CUT_TAGS

ALL_BIGWIG_CPM  = expand("{myrun}/coverage/deeptools/{sample}_CPM.bw",  sample=ALL_SAMPLES, myrun=RUNID)

########################################
# Final target list for rule all
########################################
TARGETS = []
TARGETS.extend(ALL_BIGWIG_CPM)

rule all:
    input:
        TARGETS


########################################
# Coverage: deepTools CPM
########################################
rule coverage_deeptools_cpm:
    input:
        dedup  = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bw    = "{myrun}/coverage/deeptools/{sample}_CPM.bw"
    params:
        bin_size = config['binsize'],
        smooth   = config['smooth_length']
    threads: config['THREADS']
    resources:
        mem_mb = 64000
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        bamCoverage \
          -b {input.dedup} \
          --outFileName {output.bw} \
          --outFileFormat bigwig \
          --normalizeUsing CPM \
          --binSize {params.bin_size} \
          --smoothLength {params.smooth} \
          --numberOfProcessors {threads} \
          --exactScaling
        """
