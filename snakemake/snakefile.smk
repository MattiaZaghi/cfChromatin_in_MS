#config
#configfile: "./snakemake/config_Cut_Tag.yaml"



import json

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
print(MARK_SAMPLES)
# which sample_type is used as control for calling peaks: e.g. Input, IgG...
CUT_TAG = config["c_t"]
CHIP = config["chip"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
CHIPS_SE = [sample for sample in MARK_SAMPLES if CHIP in sample]
RUNID = config["RUN_ID"]
## list BAM files
ALL_SAMPLES =  CHIPS + CUT_TAGS

BAM=expand("{myrun}/filter/samtools/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/filter/samtools/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
BED=expand("{myrun}/bed/bedtools/{sample}.bed", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG= expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
SIZE=expand("{myrun}/filter/samtools/{sample}_insert.pdf", sample = ALL_SAMPLES,myrun=RUNID)



TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(SIZE)
TARGETS.extend(BED)







ruleorder: trimming_trimmomatic > aligning_bwa > filtered_sorted_samtools > dedup_picard > filter_chr_samtools > filter_stat > coverage_deeptools > insertsize_picard > bam_to_bed 




rule all:
    input: TARGETS

import re

def get_R1(wildcards):
    sample = wildcards.sample
    files = FILES[sample.split('_')[0]][sample.split('_')[1]][sample.split('_')[2]]
    for file in files:
        if re.search(r'_R1_', file):
            return file
    raise ValueError("R1 file not found")

def get_R2(wildcards):
    sample = wildcards.sample
    files = FILES[sample.split('_')[0]][sample.split('_')[1]][sample.split('_')[2]]
    for file in files:
        if re.search(r'_R2_', file):
            return file
    raise ValueError("R2 file not found")

rule trimming_trimmomatic:
    input:
        R1 = get_R1,
        R2 = get_R2,
        adapters = config['adapters'] #"adapters/trimmomatic/adapters-pe.fa"
    output:
        Paired1       = temp("{myrun}/trimmed/trimmomatic/{sample}_R1_paired.fastq"),
        Paired2       = temp("{myrun}/trimmed/trimmomatic/{sample}_R2_paired.fastq"),
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
    """
    Align reads with bowtie2
    """
    input:
        Paired1= "{myrun}/trimmed/trimmomatic/{sample}_R1_paired.fastq",
        Paired2 = "{myrun}/trimmed/trimmomatic/{sample}_R2_paired.fastq",
    log:
        error    = "{myrun}/mapped/bowtie2/{sample}.log"
    output:
        sam = temp("{myrun}/mapped/bowtie2/{sample}.sam")
    params:
        index    = config['index_bt2_hg'] ,
        dir      = "{myrun}/mapped/bowtie2"
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
    """
    Align paired-end reads with BWA-MEM v0.7.17.
    Output: unsorted SAM (same as the old Bowtie2 rule)
    """
    input:
        Paired1 = "{myrun}/trimmed/trimmomatic/{sample}_R1_paired.fastq",
        Paired2 = "{myrun}/trimmed/trimmomatic/{sample}_R2_paired.fastq"
    output:
        sam = temp("{myrun}/mapped/bwa/{sample}.sam")
    log:
        "{myrun}/mapped/bwa/{sample}.log"
    params:
        index = config["index_bwa_hg"],    # prefix used when you ran `bwa index`
        dir   = "{myrun}/mapped/bwa"
    threads: config["THREADS"]
    resources:
        mem_mb = 64000
    conda:
        "/home/mattia/miniconda3/envs/bwa.yml"      # contains bwa 0.7.17
    shell:
        """
        mkdir -p {params.dir}
        bwa mem -M -t {threads} {params.index} {input.Paired1} {input.Paired2} > {output.sam} 2> {log}
        """


rule filtered_sorted_samtools:
    """
    1. Keep only primary alignments with MAPQ >= 10   (unique reads)
    2. Sort by coordinate
    3. Remove PCR/optical duplicates
    4. Index the BAM
    """
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
        "/home/mattia/miniconda3/envs/samtools.yml"   # samtools >=1.10 for markdup
    shell:
        """
        mkdir -p {params.dir}

        samtools view -b -F 0x904 -q 10 {input.sam} -o {output.bam}

        samtools sort -@ {threads} {output.bam} -o {output.sorted}

        samtools index -@ {threads} {output.sorted}
        """


rule dedup_picard:
    input:
        bam =  "{myrun}/sorted/samtools/{sample}.bam"
    output:
        RG = temp("{myrun}/dedup/picard/{sample}_RG.bam"),
        dedup = temp("{myrun}/dedup/picard/{sample}.bam"),
        metrics = "{myrun}/dedup/picard/{sample}.bam_metrics.txt"
    params:
        dir="{myrun}/dedup/picard",
        tmp="{myrun}/dedup/picard/tmp",
        mem="200g"
    resources:
        mem_mb=300000
    threads: config['THREADS'] 
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        picard AddOrReplaceReadGroups I={input.bam} O={output.RG} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1
        
        java -Xmx{params.mem} -jar picard MarkDuplicates I={output.RG} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp} 

        samtools index {output.dedup} -@ {threads}

        """

#  Non canonical chromosomes, Y and X unmapped and mitochondrial are removed
rule filter_chr_samtools:
    """
    Filter out Non-canonical chromosomes, Y and X  and mitochondrial
    """
    input:
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
    output:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    params:
        dir    = "{myrun}/filter/samtools"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        samtools view  -@ {threads} -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {output.filter} 
        
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
#rule down_sample:
    #input:  
        #filter ="{myrun}/filter/samtools/{sample}.bam",
        #stat ="{myrun}/filter/samtools/{sample}.flagstat"
    #output: 
        #downsample_bam="{myrun}/downsample/sambamba/{sample}.bam", 
        #downsample_bai="{myrun}/downsample/sambamba/{sample}.bam.bai"
    #resources:
        #mem_mb=64000, cpus=20
    #threads: config['THREADS']
    #run:
        #import re
        #import subprocess
        #with open (input[1], "r") as f:
            #fifth line contains the number of mapped reads
            #line = f.readlines()[4]
            #match_number = re.match(r'(\d.+) \+.+', line)
            #total_reads = int(match_number.group(1))

        #target_reads = config["target_reads"] # 15million reads  by default, set up in the config.yaml file
        #if total_reads > int(target_reads):
            #down_rate = int(target_reads)/total_reads
       # else:
            #down_rate = 1

        #shell("/home/mattia/miniconda3/envs/sambamba/bin/sambamba view -f bam -t {threads} --subsampling-seed=3 -s {rate} {inbam} |  /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} ".format(rate = down_rate, inbam = input[0], outbam = output[0]))

        #shell("/home/mattia/miniconda3/envs/samtools/bin/samtools index {outbam}".format(outbam = output[0]))

#rule downsample_stat_samtools:
    #input:
        #downsample_bam = "{myrun}/downsample/sambamba/{sample}.bam"
    #output:
        #downsample_flagstat = "{myrun}/downsample/sambamba/{sample}.flagstat"
    #resources: mem_mb=64000
    #threads: config['THREADS']
    #conda:
        #"/home/mattia/miniconda3/envs/samtools.yml"
    #shell:
        #"""        
        #samtools flagstat -@ {threads} {input.downsample_bam} > {output.downsample_flagstat} 
        
       # """

rule coverage_deeptools:
    input: 
        dedup  = "{myrun}/filter/samtools/{sample}.bam"
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

#rule R_script:
    #input: 
        #bed  = "{myrun}/dedup/picard/{sample}.bam",
        #rdata  = "Analysis/Samples/H3K27ac_ref/file.RData"  # Add your RData file here
    #output:
        #bw="{myrun}/coverage/deeptools/R/{sample}.bw"
    #params:
        #Output="Analysis/Output/H3K27ac_ref"
    #resources:
        #mem_mb=64000
    #threads: config['THREADS']
    #log: "{myrun}/coverage/deeptools/{sample}.log"
    #conda:
        #"/home/mattia/miniconda3/envs/deeptools.yml"
    #shell:
        #"""
        #Rscript --vanilla Analysis/cfChIP-seq/ProcessBEDFiles.R -r Analysis -m ${m} -BCN ${sample}
        #"""

#rule coverage_deeptools_R:
   # input: 
        #dedup  = "{myrun}/dedup/picard/{sample}.bam",
        #rdata  = "Analysis/Samples/H3K27ac_ref/file.RData"  # Add your RData file here
    #output:
        #bw="{myrun}/coverage/deeptools/R/{sample}.bw"
    #params:
        #genome_size_bp  = config['genome_size_bp'],
        #mapping_qual_bw = config['binsize'],
        #smooth= config['smooth_length']
    #resources:
       # mem_mb=64000
    #threads: config['THREADS']
    #log: "{myrun}/coverage/deeptools/{sample}.log"
    #conda:
        #"/home/mattia/miniconda3/envs/deeptools.yml"
    #shell:
        #"""
        # Use Rscript to extract the QQnorm value from the RData file
        #QQnorm=$(Rscript -e "load('{input.rdata}'); print(QQnorm)")

        # Use the extracted QQnorm value as the --scaleFactor
       # bamCoverage -b {input.dedup} --outFileName {output.bw} --scaleFactor $QQnorm --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --exactScaling
        #"""
