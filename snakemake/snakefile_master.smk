# Master Snakefile for cfChromatin Pipeline with Optional SNAPIE-Inspired Features
# This file can include all modules optionally - keep the core pipeline intact

# Configuration
configfile: "config_Cut_Tag.yaml"

import json
import os

FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())

# List all samples by sample_name, sample_type, and assay
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type+ "_" + assay)

print("Samples identified:", MARK_SAMPLES)

# which sample_type is used as control for calling peaks
CUT_TAG = config["c_t"]
CHIP = config["chip"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
ALL_SAMPLES =  CHIPS + CUT_TAGS

RUNID = config["RUN_ID"]

# ===== MAIN TARGETS (Original cfChromatin pipeline) =====
BAM = expand("{myrun}/filter/samtools/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/filter/samtools/{sample}.flagstat", sample=ALL_SAMPLES, myrun=RUNID)
BED = expand("{myrun}/bed/bedtools/{sample}.bed", sample=ALL_SAMPLES, myrun=RUNID)
ALL_BIGWIG = expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample=ALL_SAMPLES, myrun=RUNID)
SIZE = expand("{myrun}/filter/samtools/{sample}_insert.pdf", sample=ALL_SAMPLES, myrun=RUNID)
PEAKS_NARROW = expand("{myrun}/peaks/macs3/{sample}_peaks.narrowPeak", sample=ALL_SAMPLES, myrun=RUNID)
PEAKS_SUMMITS = expand("{myrun}/peaks/macs3/{sample}_summits.bed", sample=ALL_SAMPLES, myrun=RUNID)

# ===== OPTIONAL TARGETS (Extended features) =====
# Uncomment/configure these to enable additional analyses

# FastQC targets
FASTQC_ENABLED = config.get('fastqc_enabled', False)
FASTQC_TARGETS = expand("{myrun}/fastqc/{sample}", sample=ALL_SAMPLES, myrun=RUNID) if FASTQC_ENABLED else []

# Preseq (library complexity) targets
PRESEQ_TARGETS = expand("{myrun}/preseq/{sample}.lc_extrap.txt", sample=ALL_SAMPLES, myrun=RUNID) if config.get('preseq_enabled', False) else []

# Fragment length targets
FRAG_LENGTH_TARGETS = expand("{myrun}/fragments/{sample}_fraglength.txt", sample=ALL_SAMPLES, myrun=RUNID) if config.get('frag_length_enabled', False) else []

# Peak annotation targets (requires GTF)
PEAK_ANNO_TARGETS = expand("{myrun}/peaks/annotated/{sample}_peaks.annotated.txt", sample=ALL_SAMPLES, myrun=RUNID) if (config.get('peaks_annotation_enabled', False) and config.get('gtf_annotation', '')) else []

# Motif/GC targets (requires genome FASTA)
NMER = config.get('nmer_motif', 3)
MOTIF_TARGETS = expand("{myrun}/motifs/{sample}_{nmer}bp_motif.bed", sample=ALL_SAMPLES, myrun=RUNID, nmer=NMER) if (config.get('motif_gc_enabled', False) and config.get('genome_fasta', '')) else []

# Enrichment targets
ENRICH_MARKS = config.get('enrichment_marks', ['H3K4me3', 'H3K27ac'])
ENRICH_TARGETS = expand("{myrun}/enrichment/{sample}_{mark}_enrichment.txt", sample=ALL_SAMPLES, myrun=RUNID, mark=ENRICH_MARKS) if config.get('enrichment_enabled', False) else []

# SNP fingerprinting targets
SNP_TARGETS = "{myrun}/snp_fingerprint/pval_out.txt" if config.get('snp_fingerprint_enabled', False) else []

# QC report targets
QC_TARGETS = [
    expand("{myrun}/qc_reports/quality_summary.csv", myrun=RUNID),
    expand("{myrun}/qc_reports/multiqc_report.html", myrun=RUNID)
] if config.get('generate_qc_reports', False) else []

# ===== COMBINE ALL TARGETS =====
TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(SIZE)
TARGETS.extend(BED)
TARGETS.extend(PEAKS_NARROW)
TARGETS.extend(PEAKS_SUMMITS)

# Add optional targets
TARGETS.extend(FASTQC_TARGETS)
TARGETS.extend(PRESEQ_TARGETS)
TARGETS.extend(FRAG_LENGTH_TARGETS)
TARGETS.extend(PEAK_ANNO_TARGETS)
TARGETS.extend(MOTIF_TARGETS)
TARGETS.extend(ENRICH_TARGETS)
if SNP_TARGETS:
    TARGETS.append(SNP_TARGETS)
TARGETS.extend(QC_TARGETS)

ruleorder: trimming_trimmomatic > aligning_bwa > filtered_sorted_samtools > dedup_picard > filter_chr_samtools > filter_stat > coverage_deeptools > insertsize_picard > bam_to_bed > bam_to_bed_bedtools > macs3

rule all:
    input: TARGETS

import re

def get_R1(wildcards):
    sample = wildcards.sample
    files = FILES[sample.split('_')[0]][sample.split('_')[1]][sample.split('_')[2]]
    for file in files:
        if re.search(r'_R1_', file):
            return file
    raise ValueError("R1 file not found for sample: " + sample)

def get_R2(wildcards):
    sample = wildcards.sample
    files = FILES[sample.split('_')[0]][sample.split('_')[1]][sample.split('_')[2]]
    for file in files:
        if re.search(r'_R2_', file):
            return file
    raise ValueError("R2 file not found for sample: " + sample)

# ===== CORE PIPELINE RULES (Unchanged) =====

rule trimming_trimmomatic:
    input:
        R1 = get_R1,
        R2 = get_R2,
        adapters = config['adapters']
    output:
        Paired1 = temp("{myrun}/trimmed/trimmomatic/{sample}_R1_paired.fastq"),
        Paired2 = temp("{myrun}/trimmed/trimmomatic/{sample}_R2_paired.fastq"),
        Unpaired1 = temp("{myrun}/trimmed/trimmomatic/{sample}_R1_unpaired.fastq"),
        Unpaired2 = temp("{myrun}/trimmed/trimmomatic/{sample}_R2_unpaired.fastq")
    params:
        dir = "{myrun}/trimmed/trimmomatic/"
    resources:
        mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/trimmomatic.yml"
    shell:
        """
        mkdir -p {params.dir}
        trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} {output.Paired1} {output.Unpaired1} {output.Paired2} {output.Unpaired2} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule aligning_bowtie2:
    input:
        Paired1= "{myrun}/trimmed/trimmomatic/{sample}_R1_paired.fastq",
        Paired2 = "{myrun}/trimmed/trimmomatic/{sample}_R2_paired.fastq",
    log:
        error = "{myrun}/mapped/bowtie2/{sample}.log"
    output:
        sam = temp("{myrun}/mapped/bowtie2/{sample}.sam")
    params:
        index = config['index_bt2_hg'],
        dir = "{myrun}/mapped/bowtie2"
    threads: config['THREADS']
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/bowtie2.yml"
    shell:
        """
        mkdir -p {params.dir}
        bowtie2 -x {params.index} -1 {input.Paired1} -2 {input.Paired2} -S {output.sam} -p {threads} --no-mixed --no-discordant
        """

rule aligning_bwa:
    input:
        Paired1 = "{myrun}/trimmed/trimmomatic/{sample}_R1_paired.fastq",
        Paired2 = "{myrun}/trimmed/trimmomatic/{sample}_R2_paired.fastq"
    output:
        sam = temp("{myrun}/mapped/bwa/{sample}.sam")
    log:
        "{myrun}/mapped/bwa/{sample}.log"
    params:
        index = config["index_bwa_hg"],
        dir = "{myrun}/mapped/bwa"
    threads: config["THREADS"]
    resources:
        mem_mb = 64000
    conda:
        "/home/mattia/miniconda3/envs/bwa.yml"
    shell:
        """
        mkdir -p {params.dir}
        bwa mem -M -t {threads} {params.index} {input.Paired1} {input.Paired2} > {output.sam} 2> {log}
        """

rule filtered_sorted_samtools:
    input:
        sam = "{myrun}/mapped/bwa/{sample}.sam"
    output:
        bam = temp("{myrun}/view/samtools/{sample}.bam"),
        sorted = temp("{myrun}/sorted/samtools/{sample}.bam")
    params:
        dir = "{myrun}/"
    threads: config["THREADS"]
    resources:
        mem_mb = 64000
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        samtools view -b -F 0x904 -q 10 {input.sam} -o {output.bam}
        samtools sort -@ {threads} {output.bam} -o {output.sorted}
        samtools index -@ {threads} {output.sorted}
        """

rule dedup_picard:
    input:
        bam = "{myrun}/sorted/samtools/{sample}.bam"
    output:
        RG = temp("{myrun}/dedup/picard/{sample}_RG.bam"),
        dedup = temp("{myrun}/dedup/picard/{sample}.bam"),
        metrics = "{myrun}/dedup/picard/{sample}.bam_metrics.txt"
    params:
        dir="{myrun}/dedup/picard",
        tmp="{myrun}/dedup/picard/tmp"
    resources:
        mem_mb=300000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        picard AddOrReplaceReadGroups I={input.bam} O={output.RG} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1
        picard MarkDuplicates I={output.RG} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp}
        samtools index {output.dedup} -@ {threads}
        """

rule filter_chr_samtools:
    input:
        dedup = "{myrun}/dedup/picard/{sample}.bam"
    output:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    params:
        dir = "{myrun}/filter/samtools"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        samtools view -@ {threads} -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {output.filter}
        samtools index {output.filter} -@ {threads}
        """

rule filter_stat:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam",
    output:
        flagstat= "{myrun}/filter/samtools/{sample}.flagstat"
    resources: mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        samtools flagstat -@ {threads} {input.filter} > {output.flagstat}
        """

rule coverage_deeptools:
    input:
        dedup = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bw="{myrun}/coverage/deeptools/{sample}_RPKM.bw"
    params:
        genome_size_bp = config['genome_size_bp'],
        mapping_qual_bw = config['binsize'],
        norm= config['norm_method'],
        smooth= config['smooth_length']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/coverage/deeptools/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        bamCoverage -b {input.dedup} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --exactScaling
        """

rule insertsize_picard:
    input:
        dedup = "{myrun}/filter/samtools/{sample}.bam"
    output:
        metrics="{myrun}/filter/samtools/{sample}_insert.txt",
        pdf="{myrun}/filter/samtools/{sample}_insert.pdf"
    params:
        dir="{myrun}/filter/samtools/insert",
        tmp="{myrun}/filter/samtools/tmp"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        picard CollectInsertSizeMetrics I={input.dedup} O={output.metrics} H={output.pdf} M=0.5
        """

rule bam_to_bed:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bed = "{myrun}/bed/macs3/{sample}.bed"
    params:
        dir = "{myrun}/bed/macs3/"
    resources:
        mem_mb = 64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/macs3.yml"
    shell:
        """
        mkdir -p {params.dir}
        macs3 randsample -i {input.filter} -f BAMPE -p 100 -o {output.bed}
        """

rule bam_to_bed_bedtools:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bed = "{myrun}/bed/bedtools/{sample}.bed"
    params:
        dir = "{myrun}/bed/bedtools/"
    resources:
        mem_mb = 64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/macs3.yml"
    shell:
        """
        mkdir -p {params.dir}
        bedtools bamtobed -i {input.filter} -bedpe > {output.bed}
        """

rule macs3:
    input:
        treatment = "{myrun}/bed/macs3/{sample}.bed"
    output:
        peaks_narrow = "{myrun}/peaks/macs3/{sample}_peaks.narrowPeak",
        summits = "{myrun}/peaks/macs3/{sample}_summits.bed",
        peaks_xls = "{myrun}/peaks/macs3/{sample}_peaks.xls"
    params:
        outdir = "{myrun}/peaks/macs3/",
        gsize = config['genome_size_bp'],
        qvalue = config['peaks_qvalue']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/peaks/macs3/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/macs3.yml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.treatment} --name {wildcards.sample} --outdir {params.outdir} -f BEDPE --gsize {params.gsize} -q {params.qvalue} --keep-dup all 2> {log}
        """

# ===== OPTIONAL EXTENDED RULES =====
# Uncomment or include these conditionally based on config

# Include optional modules via includes
# Uncomment the modules you want to use:

# include: "fastqc.smk"  # Enable with fastqc_enabled: true
# include: "filter_properly_paired.smk"  # Enable with enable_properly_paired_filter: true
# include: "lib_complex_preseq.smk"  # Enable with preseq_enabled: true
# include: "frag_length_distribution.smk"  # Enable with frag_length_enabled: true
# include: "peaks_annotation.smk"  # Enable with peaks_annotation_enabled: true
# include: "motif_gc_content.smk"  # Enable with motif_gc_enabled: true
# include: "enrichment_analysis.smk"  # Enable with enrichment_enabled: true
# include: "snp_fingerprint.smk"  # Enable with snp_fingerprint_enabled: true
# include: "quality_reports.smk"  # Enable with generate_qc_reports: true
