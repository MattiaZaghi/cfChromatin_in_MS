#config
#configfile: "./snakemake/config_Cut_Tag.yaml"


FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())

# List all samples by sample_name, sample_type, and assay
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type+ "_" + assay)

# which sample_type is used as control for calling peaks: e.g. Input, IgG...
#CUT_TAG = config["c_t"]
#CHIP = config["chip"]
#CHIP_SE= config["chip-se"]
#CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
#CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
#CHIPS_SE = [sample for sample in MARK_SAMPLES if CHIP in sample]
#ACETYLATION=config['Acetylation']
#ACETYLATIONS = [sample for sample in MARK_SAMPLES if ACETYLATION in sample]
METHYLATION=config['Methylation']
METHYLATIONS = [sample for sample in MARK_SAMPLES if METHYLATION in sample]
#POLYCOMB=config['Polycomb']
#POLYCOMBS = [sample for sample in MARK_SAMPLES if POLYCOMB in sample]
RUNID = config["RUN_ID"]

ALL_SAMPLES = METHYLATIONS 
MULTIBAM = expand("{myrun}/coverage/deeptools/summary/multiBigwigSummary_36.npz", myrun=RUNID,sample = ALL_SAMPLES)
PCA = expand("{myrun}/coverage/deeptools/correlation/pca_36.png", sample = ALL_SAMPLES,myrun=RUNID)
CORRELATION = expand("{myrun}/coverage/deeptools/correlation/pearson_36.png", sample = ALL_SAMPLES,myrun=RUNID)


TARGETS = []
TARGETS.extend(MULTIBAM)
TARGETS.extend(PCA)
TARGETS.extend(CORRELATION)

ruleorder: multibamsummary > plot_correlation > plot_pca


rule all:
    input: TARGETS


rule multibamsummary:
    input: 
        bam =expand("{myrun}/coverage/deeptools/dedup/CPM/{sample}_CPM.bw", sample = ALL_SAMPLES,myrun=RUNID)
    output:
        multibamsummary= "{myrun}/coverage/deeptools/summary/multiBigwigSummary_36.npz",
        tab="{myrun}/coverage/deeptools/summary/multiBigwigSummary_36.tab"
    params:
        bins=config['binsize']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: 
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        multiBigwigSummary bins -b {input.bam} -o {output.multibamsummary} -bs {params.bins} -p {threads} --outRawCounts {output.tab}

        """
rule plot_correlation:
    input: 
        multibamsummary= "{myrun}/coverage/deeptools/summary/multiBigwigSummary_36.npz"
    output:
        correlation= "{myrun}/coverage/deeptools/correlation/pearson_36.png",
        cor_tab="{myrun}/coverage/deeptools/correlation/pearson_36.tab"
    params:
        plot= config['corplot'],
        stat=config['method']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: 
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        plotCorrelation --corData {input.multibamsummary} -o {output.correlation} --whatToPlot {params.plot} --corMethod {params.stat} --outFileCorMatrix {output.cor_tab}

        """


rule plot_pca:
    input: 
        multibamsummary= "{myrun}/coverage/deeptools/summary/multiBigwigSummary_36.npz"
    output:
        pca= "{myrun}/coverage/deeptools/correlation/pca_36.png",
        pca_tab="{myrun}/coverage/deeptools/correlation/pca_36.tab"
    params:
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: 
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        plotPCA --corData {input.multibamsummary} -o {output.pca} --outFileNameData {output.pca_tab}

        """