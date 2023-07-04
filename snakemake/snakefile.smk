#config
configfile: "cfChromatin_in_MS/snakemake/Cut&Tag_bulk/config.yaml"


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
RUNID = config["RUN_ID"]
print(CUT_TAGS)
print(CHIPS)
## list BAM files
ALL_SAMPLES = CUT_TAGS + CHIPS
print(ALL_SAMPLES)
BAM=expand("{myrun}/filter/samtools/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
print(BAM)
ALL_FLAGSTAT = expand("{myrun}/filter/samtools/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
print(ALL_FLAGSTAT)
ALL_DOWNSAMPLE_BAM = expand("{myrun}/downsample/sambamba/{sample}.bam", sample = ALL_SAMPLES,myrun=RUNID)
print(ALL_DOWNSAMPLE_BAM)
ALL_BIGWIG = expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
print(ALL_BIGWIG)
GOPEAKS = expand("{myrun}/peaks/gopeaks/{sample}_peaks.bed", sample = CUT_TAGS,myrun=RUNID)
print(GOPEAKS)
MACS2 = expand("{myrun}/peaks/macs2/{sample}_peaks.narrowPeak", sample = CHIPS,myrun=RUNID)
print(MACS2)

TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(ALL_DOWNSAMPLE_BAM)
TARGETS.extend(GOPEAKS)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(MACS2)
TARGETS.extend(ALL_FLAGSTAT)





ruleorder: merge_fastq > fastq_trimming > bam__bowtie2 > bam__sorted > bam__dedup > bam__filter > stat > down_sample > downsample_stat > coverage > Gopeaks > macs2

fastq_paths = FILES[sample][sample_type][assay]

rule all:
    input: TARGETS

rule merge_fastq:
    input:
        fq_1= [path for path in fastq_paths if "_R1_" in path],
        fq_2= [path for path in fastq_paths if "_R2_" in path]
    output:	
        R1 = temp("{myrun}/cat/{sample}_R1.fastq"),
        R2 = temp("{myrun}/cat/{sample}_R2.fastq")
    threads: 20
    shell:
        """
        gunzip -c {input.fq_1} > {output.R1} 
        gunzip -c {input.fq_2} > {output.R2} 
        """


rule fastq_trimming:  
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
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/trimmomatic.yml"
    shell:
        """
        mkdir -p {params.dir} 
        
        trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} {output.Paired1} {output.Unpaired1} {output.Paired2} {output.Unpaired2} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """    

rule bam__bowtie2:
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
    threads: 20
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/bowtie2.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        bowtie2 --very-sensitive-local -x {params.index} -1 {input.Paired1} -2 {input.Paired2} -S {output.sam} -p {threads} 
        """
        

rule bam__sorted:
    input:
        sam = "{myrun}/mapped/bowtie2/{sample}.sam"
    output:
        bam = temp("{myrun}/sorted/samtools/{sample}.bam")
    params:
        dir = "{myrun}/sorted/samtools/"
    threads: 20
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
        dedup = temp("{myrun}/dedup/picard/{sample}.bam"),
        bai   = temp("{myrun}/dedup/picard/{sample}.bam.bai"),
        metrics = temp("{myrun}/dedup/picard/{sample}.bam_metrics.txt")
    params:
        dir="{myrun}/dedup/picard/",
        tmp="{myrun}/dedup/picard/tmp"
    resources:
        mem_mb=140000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        java -Xmx120g -jar /home/mattia/picard.jar MarkDuplicates I={input.markdup} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp}
        
        samtools index {output.dedup} -@ {threads}
        
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
        bai="{myrun}/filter/samtools/{sample}.bam.bai"
    params:
        dir    = "{myrun}/filter/samtools/"
    resources:
        mem_mb=140000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        samtools view -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > {output.filter} -@ {threads}
        
        samtools index {output.filter} -@ {threads}
        
        """
rule stat:
    """
    Filter out Non-canonical chromosomes, Y and X  and mitochondrial
    """
    input:
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        flagstat = "{myrun}/filter/samtools/{sample}.flagstat"
    resources:
        mem_mb=140000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/samtools"
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
    threads: 20
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

        shell("/home/mattia/miniconda3/envs/sambamba view -f bam -t {threads} --subsampling-seed=3 -s {rate} {inbam} |  /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} ".format(rate = down_rate, inbam = input[0], outbam = output[0]))

        shell("/home/mattia/miniconda3/envs/samtools index {outbam}".format(outbam = output[0]))

rule downsample_stat:
    input:
        downsample_bam = "{myrun}/downsample/sambamba/{sample}.bam"
    output:
        downsample_flagstat = "{myrun}/downsample/sambamba/{sample}.flagstat"
    resources: mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/samtools"
    shell:
        """        
        /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools flagstat -@ {threads} {input.downsample_bam} > {output.downsample_flagstat} 
        
        """

rule coverage:
    input: 
        filter ="{myrun}/filter/samtools/{sample}.bam",
        bai="{myrun}/filter/samtools/{sample}.bam.bai"
    output:
        bw="{myrun}/coverage/deeptools/{sample}_RPKM.bw"
    params:
        genome_size_bp  = config['genome_size_bp'],
        mapping_qual_bw = config['binsize'],
        norm= config['norm_method'],
        smooth= config['smooth_length']
    resources:
        mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/minconda3/envs/deeptools"
    shell:
        """
        bamCoverage -b {input.filter} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.mapping_qual_bw} --smoothLength {params.smooth} --numberOfProcessors {threads} --effectiveGenomeSize {params.genome_size_bp}  --ignoreDuplicates  --skipNAs --exactScaling  
          
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
    threads: 20
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
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/macs2.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        macs2 callpeak -t {input.treatment} --name {wildcards.treatment} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --call-summits
        
        """
