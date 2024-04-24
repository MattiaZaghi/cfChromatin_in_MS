#config
#configfile: "./snakemake/config_tagAlign.yaml"


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
CHIP_SE= config["chip-se"]

CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
CHIPS_SE = [sample for sample in MARK_SAMPLES if CHIP in sample]
RUNID = config["RUN_ID"]
## list BAM files
ALL_SAMPLES = CHIPS_SE
BAM=expand("{myrun}/dedup/picard/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/dedup/picard/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG= expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BED=expand("{myrun}/bed/bedtools/{sample}.bed", sample = ALL_SAMPLES,myrun=RUNID)
ALL_SIZE=expand("{myrun}/dedup/picard/insert/{sample}_insert.pdf", sample = ALL_SAMPLES,myrun=RUNID)



TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(ALL_BED)
TARGETS.extend(ALL_SIZE)


ruleorder:  merge_fastq_SE >  fastq_trimming_SE > bam__bowtie2_SE > bam__sorted_SE > bam__dedup > coverage > stat > insertsize_picard > bam_to_bed




rule all:
    input: TARGETS


rule merge_fastq_SE:
    input:
        fq_1=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][0]
    output:
        R1=temp("{myrun}/cat/{sample}.fastq")
    log:
        R1="{myrun}/cat/{sample}.log"
    threads: config['THREADS']
    shell:
        """
        gunzip -c {input.fq_1} > {output.R1}

        """

rule fastq_trimming_SE:  
    input:
        R1       = "{myrun}/cat/{sample}.fastq", 
        adapters = config['adapters'] #"adapters/trimmomatic/adapters-pe.fa"
    output:
        Paired1       = temp("{myrun}/trimmed/trimmomatic/{sample}_trimmed.fastq")
    log: "{myrun}/trimmed/trimmomatic/{sample}.log"
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
        
        trimmomatic SE -threads {threads} -phred33 {input.R1} {output.Paired1} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

        """    

rule bam__bowtie2_SE:
    """
    Align reads with bowtie2
    """
    input:
        Paired1= "{myrun}/trimmed/trimmomatic/{sample}_trimmed.fastq"
    output:
        sam = temp("{myrun}/mapped/bowtie2/{sample}.sam")
    params:
        index    = config['index_bt2_hg'] ,
        dir      = "{myrun}/mapped/bowtie2/"
    log: "{myrun}/mapped/bowtie2/SE/{sample}.log"
    threads: config['THREADS']
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/bowtie2.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        bowtie2 -x {params.index} -U {input.Paired1} -S {output.sam} -p {threads} --no-mixed --no-discordant

        """
rule bam__sorted_SE:
    input:
        sam = "{myrun}/mapped/bowtie2/{sample}.sam"
    output:
        bam = "{myrun}/sorted/samtools/{sample}.bam"
    params:
        dir = "{myrun}/sorted/samtools/"
    threads: config['THREADS'] 
    log: "{myrun}/sorted/samtools/{sample}.log"
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        samtools sort -o {output.bam} -O bam {input.sam} -@ {threads} 
        
        samtools index {output.bam} -@ {threads}

        """


rule bam__dedup:
    input:
        markdup = "{myrun}/sorted/samtools/{sample}.bam"
    output:
        dedup = "{myrun}/dedup/picard/{sample}.bam",
        bai   = "{myrun}/dedup/picard/{sample}.bam.bai",
        metrics = "{myrun}/dedup/picard/{sample}.bam_metrics.txt"
    params:
        dir="{myrun}/dedup/picard/",
        tmp="{myrun}/dedup/picard/tmp"
    resources:
        mem_mb=140000
    threads: config['THREADS'] 
    log: "{myrun}/dedup/picard/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        picard MarkDuplicates I={input.markdup} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp} 

        samtools index {output.dedup} -@ {threads}

        """
rule coverage:
    input: 
        dedup = "{myrun}/dedup/picard/{sample}.bam"
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
    log: "{myrun}/coverage/deeptools/dedup/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        bamCoverage -b {input.dedup} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads}  --exactScaling  

        """
#  Non canonical chromosomes, Y and X unmapped and mitochondrial are removed
rule bam__filter:
    """
    Filter out Non-canonical chromosomes, Y and X  and mitochondrial
    """
    input:
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
    output:
        filter = "{myrun}/filter/samtools/{sample}.bam",
        bai="{myrun}/filter/samtools/{sample}.bam.bai",
        #insert_size_metrics="{myrun}/filter/samtools/{sample}_insert.txt",
        #insert_size_histogram="{myrun}/filter/samtools/{sample}_insert.pdf"
    params:
        dir    = "{myrun}/filter/samtools/"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    log: "{myrun}/filter/samtools/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        samtools view  -@ {threads} -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > {output.filter} 
        
        samtools index {output.filter} -@ {threads}
        
        """
rule stat:
    input:
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
    output:
        flagstat = "{myrun}/dedup/picard/{sample}.flagstat"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """        
        samtools flagstat -@ {threads} {input.dedup} > {output.flagstat} 
    
        """

rule insertsize_picard:
    input:
        dedup = "{myrun}/dedup/picard/{sample}.bam"
    output:
        metrics="{myrun}/dedup/picard/insert/{sample}_insert.txt",
        pdf="{myrun}/dedup/picard/insert/{sample}_insert.pdf"
    params:
        dir="{myrun}/dedup/picard/insert",
        tmp="{myrun}/dedup/picard/tmp"
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
        filter = "{myrun}/dedup/picard/{sample}.bam"
    output:
        bed="{myrun}/bed/bedtools/{sample}.bed"
    params:
        dir  = "{myrun}/bed/bedtools/"
    resources:
        mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/bedtools.yml"
    shell:
        """

        mkdir -p {params.dir}

        bedtools bamtobed -i {input.filter}  > {output.bed}

        """
        
        