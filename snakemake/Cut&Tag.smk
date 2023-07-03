#config
configfile: "cfChromatin_in_MS/snakemake/Cut&Tag_bulk/config.yaml"


FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        MARK_SAMPLES.append(sample + "_" + assay)


# which sample_type is used as control for calling peaks: e.g. Input, IgG...
CUT_TAG = config["c_t"]
CHIP = config["chip"]

## list BAM files
CUT_TAG = expand("filter/{sample}.bam", sample=CUT_TAG)
CHIP = expand("filter/{sample}.bam", sample=CHIP)

ALL_FASTQC  = expand("02fqc/{sample}_fastqc.zip", sample = ALL_BAM)
ALL_SAMPLES = CUT_TAG + CHIP
ALL_DOWNSAMPLE_BAM = expand("downsample/{sample}_downsample.bam", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("flagstat/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_BIGWIG = expand("coverage/{sample}_RPKM.bw", sample = ALL_SAMPLES)
GOPEAKS = expand("coverage/{sample}_peaks.bed", sample = CUT_TAG)
MACS2 = expand("coverage/{sample}_peaks.narrowPeak", sample = CHIP)

TARGETS = []
TARGETS.extend(ALL_FASTQC)
TARGETS.extend(ALL_BAM)
TARGETS.extend(ALL_DOWNSAMPLE_BAM)
TARGETS.extend(GOPEAKS)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(MACS2)
TARGETS.extend(ALL_FLAGSTAT)



RUNID = config["RUN_ID"]

ruleorder: fastq_trimming__ > bam__bowtie2 > bam__sorted > bam__markdup > bam__dedup > bam__filter > stat > bam__bigWig > peak_calling_ATAC > Gopeaks_mature > Gopeaks_P4 #> down_sample > downsample_stat > bam__bigWig_downsample > peak_calling_ATAC_downsample > Gopeaks_P4_downsample > Gopeaks_mature_downsample


rule all:
    input: TARGETS

rule merge_fastqs:
    input: get_fastq
    output: "01seq/{sample}.fastq"
    log: "00log/{sample}_unzip"
    threads: CLUSTER["merge_fastqs"]["cpu"]
    params: jobname = "{sample}"
    message: "merging fastqs gunzip -c {input} > {output}"
    shell: "gunzip -c {input} > {output} 2> {log}"

rule merge_fastqs:  
    input:
        R1       = "/date/gcb/GCB_MK/P29054/P29054_1005/02-FASTQ/230510_A00187_0957_BH2CM3DRX3/{mysample}_R1_001.fastq.gz", 
        R2       = "/date/gcb/GCB_MK/P29054/P29054_1005/02-FASTQ/230510_A00187_0957_BH2CM3DRX3/{mysample}_R2_001.fastq.gz",
        adapters = config['adapters'] #"adapters/trimmomatic/adapters-pe.fa"
    output:
        Paired1       = temp("{myrun}/trimmed/trimmomatic/{mysample}_R1_paired.fastq.gz"),
        Paired2       = temp("{myrun}/trimmed/trimmomatic/{mysample}_R2_paired.fastq.gz"),
        Unpaired1 = temp("{myrun}/trimmed/trimmomatic/{mysample}_R1_unpaired.fastq.gz"),
        Unpaired2 = temp("{myrun}/trimmed/trimmomatic/{mysample}_R2_unpaired.fastq.gz")
    log:
        main     = "{myrun}/trimmed/trimmomatic/{mysample}_trim.log",
        out      = "{myrun}/trimmed/trimmomatic/{mysample}_trimout.log"
    params:
        dir = "trimmed/trimmomatic/"
    resources:
        mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/trimmomatic.yml"
    shell:
        """
        mkdir -p {params.dir} 
        
        trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} {output.Paired1} {output.Unpaired1} {output.Paired2} {output.Unpaired2} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > "{log.out}"
        """    


rule fastq_trimming__:  
    input:
        R1       = "/date/gcb/GCB_MK/P29054/P29054_1005/02-FASTQ/230510_A00187_0957_BH2CM3DRX3/{mysample}_R1_001.fastq.gz", 
        R2       = "/date/gcb/GCB_MK/P29054/P29054_1005/02-FASTQ/230510_A00187_0957_BH2CM3DRX3/{mysample}_R2_001.fastq.gz",
        adapters = config['adapters'] #"adapters/trimmomatic/adapters-pe.fa"
    output:
        Paired1       = temp("{myrun}/trimmed/trimmomatic/{mysample}_R1_paired.fastq.gz"),
        Paired2       = temp("{myrun}/trimmed/trimmomatic/{mysample}_R2_paired.fastq.gz"),
        Unpaired1 = temp("{myrun}/trimmed/trimmomatic/{mysample}_R1_unpaired.fastq.gz"),
        Unpaired2 = temp("{myrun}/trimmed/trimmomatic/{mysample}_R2_unpaired.fastq.gz")
    log:
        main     = "{myrun}/trimmed/trimmomatic/{mysample}_trim.log",
        out      = "{myrun}/trimmed/trimmomatic/{mysample}_trimout.log"
    params:
        dir = "trimmed/trimmomatic/"
    resources:
        mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/trimmomatic.yml"
    shell:
        """
        mkdir -p {params.dir} 
        
        trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} {output.Paired1} {output.Unpaired1} {output.Paired2} {output.Unpaired2} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > "{log.out}"
        """    

rule bam__bowtie2:
    """
    Align reads with bowtie2
    """
    input:
        R1       = "{myrun}/trimmed/trimmomatic/{mysample}_R1_paired.fastq.gz",
        R2       = "{myrun}/trimmed/trimmomatic/{mysample}_R2_paired.fastq.gz",
    log:
        error    = "{myrun}/mapped/bowtie2/{mysample}.log"
    output:
        sam = temp("{myrun}/mapped/bowtie2/{mysample}.sam")
    params:
        index    = config['index_bt2_hg'] ,
        dir      = "{myrun}/mapped/bowtie2"
    resources:
        mem_mb=64000
    threads: 20
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/bowtie2.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        bowtie2 --very-sensitive-local -x {params.index} -1 {input.R1} -2 {input.R2} -S {output.sam} -p {threads} > "{log.error}"
        """
        

rule bam__sorted:
    input:
        sam = "mapped/bowtie2/{mysample}.sam"
    output:
        bam = temp("mapped/bowtie2/sorted/samtools/{mysample}.bam")
    params:
        dir = "mapped/bowtie2/sorted/samtools"
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
        markdup = "mapped/bowtie2/sorted/samtools/{mysample}.bam"
    output:
        dedup = temp("mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam"),
        bai   = temp("mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam.bai"),
        metrics = temp("mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam_metrics.txt")
    params:
        dir="mapped/bowtie2/sorted/samtools/dedup/",
        tmp = "/proj/tmp/tmp_MZ/tmp"
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
        dedup  = "mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam"
    output:
        filter = "mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam",
        bai="mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.bai"
    params:
        dir    = "mapped/bowtie2/sorted/samtools/dedup/filter"
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
        filter = "mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam"
    output:
        flagstat = "mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.flagstat"
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
        filter ="mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam",
        stat ="mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.flagstat"
    output: 
        downsample_bam="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam", 
        downsample_bai="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam.bai"
    resources:
        mem_mb=64000, cpus=20
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/sambamba.yml"
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
            #down_rate = int(target_reads)/total_reads
        else:
            down_rate = 1

        shell("sambamba view -f bam -t 36 --subsampling-seed=3 -s {rate} {inbam} |  /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} ".format(rate = down_rate, inbam = input[0], outbam = output[0]))

        shell("/home/mattia/miniconda3/envs/samtools index {outbam}".format(outbam = output[0]))

rule downsample_stat:
    input:
        downsample_bam = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam"
    output:
        downsample_flagstat = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam.flagstat"
    resources: mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/miniconda3/envs/samtools"
    shell:
        """        
        /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools flagstat -@ 20 {input.downsample_bam} > {output.downsample_flagstat} 
        
        """

rule bam__bigWig:
    input: 
        bam="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam",
        bai="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.bai"
    output:
        bw="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}_RPKM.bw"
    params:
        ref_genome_fa   =  config['ref_genome_fa'],
        # blacklist       =  config['blacklist'],
        genome_size_bp  = config['genome_size_bp'],
        mapping_qual_bw = 10
    resources:
        mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/minconda3/envs/deeptools"
    shell:
        """
        bamCoverage -b {input.bam} --outFileName {output.bw} --normalizeUsing RPKM --binSize 10 --smoothLength 300 --numberOfProcessors 10 --effectiveGenomeSize {params.genome_size_bp}  --ignoreDuplicates  --skipNAs --exactScaling  
          
        """

rule Gopeaks:
    input:
        downsample_treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{pups}.bam",
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{igp4}.bam",igp4=IGP4, myrun=RUNID)
    output:
        downsample_peaks_pups =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/Gopeaks/{pups}_peaks.bed"
    params:
        peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/Gopeaks",
        qvalue     = config['peaks_qvalue'],
        outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/Gopeaks"
    resources:
        mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/minconda3/envs/gopeaks.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        gopeaks -b {input.downsample_treatment} -o {params.outdir}/{wildcards.pups} -p {params.qvalue} 
        
        """

rule macs2:
    input:
         treatment = "mapped/bowtie2/sorted/samtools/dedup/filter/{mysample_narrow}.bam",
    output:
        peaks_narrow =  "peaks/macs2/{mysample}_peaks.narrowPeak",
        bed   =  "peaks/macs2/{mysample}_summits.bed"
    params:
        peaks_dir  = "peaks/macs2/",
        gsize      = config['genome_size_bp'],
        qvalue     = config['peaks_qvalue'],
        outdir     = "peaks/macs2/"
    resources: mem_mb=64000
    threads: 20
    conda:
        "/home/mattia/minconda3/envs/macs2.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        /proj/tmp/tmp_MZ/anaconda3/envs/macs2/bin/macs2 callpeak -t {input.treatment} --name {wildcards.mysample_narrow} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --call-summits
        
        """
