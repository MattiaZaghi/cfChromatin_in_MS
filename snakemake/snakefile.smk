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

BAM=expand("{myrun}/dedup/picard/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/dedup/picard/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
FILTER expand("{myrun}/filter/samtools/{sample}.bam", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG_DEDUP = expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
BED=expand("{myrun}/bed/bedtools/{sample}.bed", sample = ALL_SAMPLES,myrun=RUNID)
#ALL_BIGWIG= expand("{myrun}/coverage/deeptools/CPM/{sample}_CPM.bw", sample = ALL_SAMPLES,myrun=RUNID)
#GOPEAKS = expand("{myrun}/peaks/gopeaks/{sample}_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
#GOPEAKS_BROAD = expand("{myrun}/peaks/gopeaks/{sample}_broad_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
#MACS2 = expand("{myrun}/peaks/macs2/{sample}_peaks.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
#MACS2_BROAD = expand("{myrun}/peaks/macs2/{sample}_peaks.broadPeak", sample = ALL_SAMPLES,myrun=RUNID)
SIZE=expand("{myrun}/dedup/picard/insert/{sample}_insert.pdf", sample = ALL_SAMPLES,myrun=RUNID)



TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(FILTER)
#TARGETS.extend(GOPEAKS)
#TARGETS.extend(GOPEAKS_BROAD)
TARGETS.extend(ALL_BIGWIG_DEDUP)
#TARGETS.extend(ALL_BIGWIG)
#TARGETS.extend(MACS2)
#TARGETS.extend(MACS2_BROAD)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(SIZE)
TARGETS.extend(BED)






ruleorder: merge_fastq >  trimming_trimmomatic >  aligning_bowtie2 >  sorted_samtools > dedup_picard > coverage_deeptools > insertsize_picard > bam_to_bed




rule all:
    input: TARGETS

rule merge_fastq:
    input:
        fq_1_L001=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][0],
        fq_1_L002=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][2],
        fq_2_L001=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][1],
        fq_2_L002=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][3]
    output:
        R1=temp("{myrun}/cat/{sample}_R1.fastq"),
        R2=temp("{myrun}/cat/{sample}_R2.fastq")
    threads: config['THREADS']
    shell:
        """
        gunzip -c {input.fq_1_L001} {input.fq_1_L002} > {output.R1}

        gunzip -c {input.fq_2_L001} {input.fq_2_L002} > {output.R2}

        """


rule trimming_trimmomatic:  
    input:
        R1       = "{myrun}/cat/{sample}_R1.fastq", 
        R2       = "{myrun}/cat/{sample}_R2.fastq",
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


rule sorted_samtools:
    input:
        sam = "{myrun}/mapped/bowtie2/{sample}.sam"
    output:
        bam = temp("{myrun}/sorted/samtools/{sample}.bam")
    params:
        dir = "{myrun}/sorted/samtools/"
    threads: config['THREADS'] 
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

rule dedup_picard:
    input:
        bam =  "{myrun}/sorted/samtools/{sample}.bam"
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
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        picard MarkDuplicates I={input.bam} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp} 

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
        filter = "{myrun}/filter/samtools/{sample}.bam",
        bai="{myrun}/filter/samtools/{sample}.bam.bai"
    params:
        dir    = "{myrun}/filter/samtools/"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        samtools view  -@ {threads} -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM > {output.filter} 
        
        samtools index {output.filter} -@ {threads}
        
        """
rule filter_stat:
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam",
    output:
        downsample_filter= "{myrun}/filter/samtools/{sample}.flagstat"
    resources: mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """        
        samtools flagstat -@ {threads} {input.downsample_bam} > {output.downsample_flagstat} 
        
        """
rule stat_samtools:
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
rule down_sample:
    input:  
        filter ="{myrun}/filter/samtools/{sample}.bam",
        stat ="{myrun}/filter/samtools/{sample}.flagstat"
    output: 
        downsample_bam="{myrun}/downsample/sambamba/{sample}.bam", 
        downsample_bai="{myrun}/downsample/sambamba/{sample}.bam.bai"
    resources:
        mem_mb=64000, cpus=20
    threads: config['THREADS']
    run:
        import re
        import subprocess
        with open (input[1], "r") as f:
            #fifth line contains the number of mapped reads
            line = f.readlines()[4]
            match_number = re.match(r'(\d.+) \+.+', line)
            total_reads = int(match_number.group(1))

        target_reads = config["target_reads"] # 15million reads  by default, set up in the config.yaml file
        if total_reads > int(target_reads):
            down_rate = int(target_reads)/total_reads
        else:
            down_rate = 1

        shell("/home/mattia/miniconda3/envs/sambamba/bin/sambamba view -f bam -t {threads} --subsampling-seed=3 -s {rate} {inbam} |  /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} ".format(rate = down_rate, inbam = input[0], outbam = output[0]))

        shell("/home/mattia/miniconda3/envs/samtools/bin/samtools index {outbam}".format(outbam = output[0]))

rule downsample_stat_samtools:
    input:
        downsample_bam = "{myrun}/downsample/sambamba/{sample}.bam"
    output:
        downsample_flagstat = "{myrun}/downsample/sambamba/{sample}.flagstat"
    resources: mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """        
        samtools flagstat -@ {threads} {input.downsample_bam} > {output.downsample_flagstat} 
        
        """

rule coverage_deeptools:
    input: 
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
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
rule coverage_deeptools_filter:
    input: 
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
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
rule R_script:
    input: 
        bed  = "{myrun}/dedup/picard/{sample}.bam",
        rdata  = "Analysis/Samples/H3K27ac_ref/file.RData"  # Add your RData file here
    output:
        bw="{myrun}/coverage/deeptools/R/{sample}.bw"
    params:
        Output="Analysis/Output/H3K27ac_ref"
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/coverage/deeptools/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        Rscript --vanilla Analysis/cfChIP-seq/ProcessBEDFiles.R -r Analysis -m ${m} -BCN ${sample}
        """

rule coverage_deeptools_R:
    input: 
        dedup  = "{myrun}/dedup/picard/{sample}.bam",
        rdata  = "Analysis/Samples/H3K27ac_ref/file.RData"  # Add your RData file here
    output:
        bw="{myrun}/coverage/deeptools/R/{sample}.bw"
    params:
        genome_size_bp  = config['genome_size_bp'],
        mapping_qual_bw = config['binsize'],
        smooth= config['smooth_length']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/coverage/deeptools/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        # Use Rscript to extract the QQnorm value from the RData file
        QQnorm=$(Rscript -e "load('{input.rdata}'); print(QQnorm)")

        # Use the extracted QQnorm value as the --scaleFactor
        bamCoverage -b {input.dedup} --outFileName {output.bw} --scaleFactor $QQnorm --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --exactScaling
        """
