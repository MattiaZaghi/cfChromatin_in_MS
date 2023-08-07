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
CHIP_SE= config["chip-se"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
CHIPS_SE = [sample for sample in MARK_SAMPLES if CHIP in sample]
RUNID = config["RUN_ID"]
## list BAM files
ALL_SAMPLES = CHIPS
BAM=expand("{myrun}/filter/samtools/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/filter/samtools/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
ALL_DOWNSAMPLE_BAM = expand("{myrun}/downsample/sambamba/{sample}.bam", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG_SORTED = expand("{myrun}/coverage/deeptools/sorted/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG_DEDUP = expand("{myrun}/coverage/deeptools/dedup/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG= expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
GOPEAKS = expand("{myrun}/peaks/gopeaks/{sample}_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
GOPEAKS_BROAD = expand("{myrun}/peaks/gopeaks/{sample}_broad_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
MACS2 = expand("{myrun}/peaks/macs2/{sample}_peaks.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
MACS2_BROAD = expand("{myrun}/peaks/macs2/{sample}_peaks.broadPeak", sample = ALL_SAMPLES,myrun=RUNID)



TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(GOPEAKS)
TARGETS.extend(GOPEAKS_BROAD)
TARGETS.extend(ALL_BIGWIG_SORTED)
TARGETS.extend(ALL_BIGWIG_DEDUP)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(MACS2)
TARGETS.extend(MACS2_BROAD)
TARGETS.extend(ALL_FLAGSTAT)





ruleorder:  bam__bowtie2_tag > bam__sorted_tag > coverage_sorted > bam__dedup > coverage_dedup > bam__filter > stat > down_sample > downsample_stat > coverage > Gopeaks > Gopeaks_broad > macs2 > macs2_broad 




rule all:
    input: TARGETS


rule bam__bowtie2_tag:
    """
    Align reads with bowtie2
    """
    input:
        tag=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][0]
    output:
        sam = temp("{myrun}/mapped/bowtie2/{sample}.sam")
    params:
        index    = config['index_bt2_hg'] ,
        dir      = "{myrun}/mapped/bowtie2/"
    threads: config['THREADS']
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/bowtie2.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        bowtie2 --very-sensitive-local -x {params.index} -U {input.tag} -S {output.sam} -p {threads} --no-mixed --no-discordant

        """
rule bam__sorted_tag:
    input:
        sam = "{myrun}/mapped/bowtie2/{sample}.sam"
    output:
        bam = "{myrun}/sorted/samtools/{sample}.bam",
        flagstat = "{myrun}/sorted/samtools/{sample}.flagstat",
        #insert_size_metrics="{myrun}/sorted/samtools/{sample}_insert.txt",
        #insert_size_histogram="{myrun}/sorted/samtools/{sample}_insert.pdf"
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

        samtools flagstat -@ {threads} {output.bam} > {output.flagstat} 

        """

rule coverage_sorted:
    input: 
        bam = "{myrun}/sorted/samtools/{sample}.bam",
    output:
        bw="{myrun}/coverage/deeptools/sorted/{sample}_RPKM.bw"
    params:
        genome_size_bp  = config['genome_size_bp'],
        mapping_qual_bw = config['binsize'],
        norm= config['norm_method'],
        smooth= config['smooth_length']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/coverage/deeptools/sorted/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        bamCoverage -b {input.bam} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --effectiveGenomeSize {params.genome_size_bp}  --exactScaling  

        """

rule bam__dedup:
    input:
        markdup = "{myrun}/sorted/samtools/{sample}.bam"
    output:
        dedup = "{myrun}/dedup/picard/{sample}.bam",
        bai   = "{myrun}/dedup/picard/{sample}.bam.bai",
        metrics = "{myrun}/dedup/picard/{sample}.bam_metrics.txt",
        flagstat = "{myrun}/dedup/picard/{sample}.flagstat",
        #insert_size_metrics="{myrun}/dedup/picard/{sample}_insert.txt",
        #insert_size_histogram="{myrun}/dedup/picard/{sample}_insert.pdf"
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
        
        samtools flagstat -@ {threads} {output.dedup} > {output.flagstat}

        """
rule coverage_dedup:
    input: 
        dedup = "{myrun}/dedup/picard/{sample}.bam"
    output:
        bw="{myrun}/coverage/deeptools/dedup/{sample}_RPKM.bw"
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
        filter = "{myrun}/filter/samtools/{sample}.bam",
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
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

rule downsample_stat:
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
        /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools flagstat -@ {threads} {input.downsample_bam} > {output.downsample_flagstat} 
        
        """

rule coverage:
    input: 
        filter ="{myrun}/filter/samtools/{sample}.bam"
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
        treatment = "{myrun}/filter/samtools/{sample}.bam",
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
    
        gopeaks -b {input.treatment} -o {params.outdir}/{wildcards.sample} -p {params.qvalue} 

        """

rule Gopeaks_broad:
    input:
        treatment = "{myrun}/filter/samtools/{sample}.bam",
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
    
        gopeaks -b {input.treatment} -o {params.outdir}/{wildcards.sample} -p {params.qvalue} --broad 

        """

rule macs2:
    input:
         treatment = "{myrun}/filter/samtools/{sample}.bam",
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
    
        macs2 callpeak -t {input.treatment} --name {wildcards.sample} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --call-summits
        
        """

rule macs2_broad:
    input:
         treatment = "{myrun}/filter/samtools/{sample}.bam",
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
    
        macs2 callpeak -t {input.treatment} --name {wildcards.sample} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --broad 
        
        """
        
        