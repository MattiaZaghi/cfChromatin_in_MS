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
CUT_TAG = config["c_t"]
CHIP = config["chip"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
CHIPS_SE = [sample for sample in MARK_SAMPLES if CHIP in sample]
RUNID = config["RUN_ID"]
## list BAM files
ALL_SAMPLES =  CUT_TAGS

ALL_BIGWIG= expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
GOPEAKS = expand("{myrun}/peaks/gopeaks/{sample}_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
GOPEAKS_BROAD = expand("{myrun}/peaks/gopeaks/{sample}_broad_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
MACS2 = expand("{myrun}/peaks/macs2/{sample}_peaks.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
MACS2_BROAD = expand("{myrun}/peaks/macs2/{sample}_peaks.broadPeak", sample = ALL_SAMPLES,myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/filter/samtools/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)




TARGETS = []
TARGETS.extend(GOPEAKS)
TARGETS.extend(GOPEAKS_BROAD)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(MACS2)
TARGETS.extend(MACS2_BROAD)
TARGETS.extend(ALL_FLAGSTAT)





ruleorder: mv_bam > coverage > Gopeaks > Gopeaks_broad > macs2 > macs2_broad > stat


rule all:
    input: TARGETS


rule mv_bam:
    input:
        bam=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][0]
    output:
        filter = temp("{myrun}/filter/samtools/{sample}.bam")
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        cp {input.bam} {output.filter}

        samtools index {output.filter}

        """
rule coverage:
    input: 
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bw="{myrun}/coverage/deeptools/{sample}_RPKM.bw"
    params:
        genome_size_bp  = config['genome_size_bp'],
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
        bamCoverage -b {input.filter} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --exactScaling  

        """
rule Gopeaks:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{igp4}.bam",igp4=IGP4, myrun=RUNID)
    output:
        peaks =  "{myrun}/peaks/gopeaks/{sample}_peaks.bed"
    params:
        peaks_dir  = "{myrun}/peaks/gopeaks/",
        qvalue     = config['peaks_qvalue'],
        outdir     = "{myrun}/peaks/gopeaks/"
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/peaks/gopeaks/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/gopeaks.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        gopeaks -b {input.filter} -o {params.outdir}/{wildcards.sample} -p {params.qvalue} 

        """

rule Gopeaks_broad:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{igp4}.bam",igp4=IGP4, myrun=RUNID)
    output:
        peaks =  "{myrun}/peaks/gopeaks/{sample}_broad_peaks.bed"
    params:
        peaks_dir  = "{myrun}/peaks/gopeaks/",
        qvalue     = config['peaks_qvalue'],
        outdir     = "{myrun}/peaks/gopeaks/"
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/peaks/gopeaks/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/gopeaks.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        gopeaks -b {input.filter} -o {params.outdir}/{wildcards.sample} -p {params.qvalue} --broad 

        """

rule macs2:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        peaks_narrow =  "{myrun}/peaks/macs2/{sample}_peaks.narrowPeak",
        bed   =  "{myrun}/peaks/macs2/{sample}_summits.bed"
    params:
        peaks_dir  = "{myrun}/peaks/macs2/",
        gsize      = config['genome_size_bp'],
        qvalue     = config['peaks_qvalue'],
        outdir     = "{myrun}/peaks/macs2/"
    resources: mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/peaks/macs2/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/macs2.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        macs2 callpeak -t {input.filter} --name {wildcards.sample} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --call-summits
        
        """

rule macs2_broad:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        peaks_narrow =  "{myrun}/peaks/macs2/{sample}_peaks.broadPeak",
        bed   =  "{myrun}/peaks/macs2/{sample}_peaks.gappedPeak"
    params:
        peaks_dir  = "{myrun}/peaks/macs2/",
        gsize      = config['genome_size_bp'],
        qvalue     = config['peaks_qvalue'],
        outdir     = "{myrun}/peaks/macs2/"
    resources: mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/peaks/macs2/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/macs2.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        macs2 callpeak -t {input.filter} --name {wildcards.sample} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --broad 
        
        """
        
rule stat:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        flagstat = "{myrun}/filter/samtools/{sample}.flagstat"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """        
        samtools flagstat -@ {threads} {input.filter} > {output.flagstat} 
    
        """