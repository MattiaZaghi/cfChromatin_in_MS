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
ACETYLATION=config['Acetylation']
ACETYLATIONS = [sample for sample in MARK_SAMPLES if ACETYLATION in sample]
#METHYLATION=config['Methylation']
#METHYLATIONS = [sample for sample in MARK_SAMPLES if METHYLATION in sample]
#POLYCOMB=config['Polycomb']
#POLYCOMBS = [sample for sample in MARK_SAMPLES if POLYCOMB in sample]
RUNID = config["RUN_ID"]

ALL_SAMPLES = ACETYLATIONS
MULTIBAM = expand("{myrun}/deeptools/summary/multiBamSummary_TSS_K27ac.npz", myrun=RUNID,sample = ALL_SAMPLES)
PCA = expand("{myrun}/deeptools/correlation/pca_TSS_K27ac.png", sample = ALL_SAMPLES,myrun=RUNID)
CORRELATION = expand("{myrun}/deeptools/correlation/pearson_TSS_K27ac.png", sample = ALL_SAMPLES,myrun=RUNID)


TARGETS = []
TARGETS.extend(MULTIBAM)
TARGETS.extend(PCA)
TARGETS.extend(CORRELATION)


ruleorder: multibamsummary > plot_correlation > plot_pca


rule all:
    input: TARGETS


rule multibamsummary:
    input: 
        bam =expand("{myrun}/dedup/picard/{sample}.bam", sample = ALL_SAMPLES,myrun=RUNID),
        bed=expand("{myrun}/TSS_hg38_list.bed", sample = ALL_SAMPLES,myrun=RUNID),
    output:
        multibamsummary= "{myrun}/deeptools/summary/multiBamSummary_TSS_K27ac.npz",
        tab="{myrun}/deeptools/summary/multiBamSummary_TSS_K27ac.tab",
        norm_factor="{myrun}/deeptools/summary/norm_factor_TSS_K27ac.tab"
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
        multiBamSummary BED-file -b {input.bam} -o {output.multibamsummary}  -p {threads} --outRawCounts {output.tab} --BED {input.bed} --scalingFactors {output.norm_factor} 

        """
rule plot_correlation:
    input: 
        multibamsummary= "{myrun}/deeptools/summary/multiBamSummary_TSS_K27ac.npz"
    output:
        correlation= "{myrun}/deeptools/correlation/pearson_TSS_K27ac.png",
        cor_tab="{myrun}/deeptools/correlation/pearson_TSS_K27ac.tab"
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
        multibamsummary= "{myrun}/deeptools/summary/multiBamSummary_TSS_K27ac.npz"
    output:
        pca= "{myrun}/deeptools/correlation/pca_TSS_K27ac.png",
        pca_tab="{myrun}/deeptools/correlation/pca_TSS_K27ac.tab"
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