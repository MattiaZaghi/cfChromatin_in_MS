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

ALL_SAMPLES = CUT_TAG + CHIP
ALL_BAM     = CONTROL_BAM + CASE_BAM
ALL_DOWNSAMPLE_BAM = expand("downsample/{sample}-downsample.sorted.bam", sample = ALL_SAMPLES)
ALL_FASTQC  = expand("02fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("flagstat/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_BIGWIG = expand("Coverage/{sample}_RPKM.bw", sample = ALL_SAMPLES)
GOPEAKS = expand("Coverage/{sample}_peaks.bed", sample = CUT_TAG)
MACS2 = expand("Coverage/{sample}_peaks.narrowPeak", sample = CHIP)

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


rule fastq_trimming__:  
    input:
        R1       = "/date/gcb/GCB_MK/P29054/P29054_1005/02-FASTQ/230510_A00187_0957_BH2CM3DRX3/zcat/{mysample}_R1_001.fastq.gz", 
        R2       = "/date/gcb/GCB_MK/P29054/P29054_1005/02-FASTQ/230510_A00187_0957_BH2CM3DRX3/zcat/{mysample}_R2_001.fastq.gz",
        adapters = config['adapters'] #"adapters/trimmomatic/adapters-pe.fa"
    output:
        Paired1       = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_R1_paired.fastq.gz"),
        Paired2       = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_R2_paired.fastq.gz"),
        Unpaired1 = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_R1_unpaired.fastq.gz"),
        Unpaired2 = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_R2_unpaired.fastq.gz")
    log:
        main     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_trim.log",
        out      = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_trimout.log"
    params:
        dir = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/"
    resources:
        mem_mb=64000
    threads: 20
    #conda:
        #"trimmomatic"
    shell:
        """
        mkdir -p {params.dir} 
        
        /home/mattia/miniconda3/envs/trimmomatic/bin/trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} {output.Paired1} {output.Unpaired1} {output.Paired2} {output.Unpaired2} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > "{log.out}"
        """    

rule bam__bowtie2:
    """
    Align reads with bowtie2
    """
    input:
        R1       = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_R1_paired.fastq.gz",
        R2       = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/{mysample}_R2_paired.fastq.gz",
    log:
        error    = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/{mysample}.log"
    output:
        sam = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/{mysample}.sam")
    params:
        index    = "/home/mattia/Genomes/hg38/fa/hg38",
        dir      = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2"
    resources:
        mem_mb=64000
    threads: 20
    resources:
        mem_mb=64000
    #conda:
        #"bowtie2"
    shell:
        """
        mkdir -p {params.dir}
        
        /home/mattia/miniconda3/envs/bowtie2/bin/bowtie2 --very-sensitive-local -x {params.index} -1 {input.R1} -2 {input.R2} -S {output.sam} -p {threads} > "{log.error}"
        """
        

rule bam__sorted:
    input:
        sam = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/{mysample}.sam"
    output:
        bam = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/{mysample}.bam")
    params:
        dir = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools"
    threads: 20
    resources:
        mem_mb=64000
    #conda:
        #"samtools"
    shell:
        """
        mkdir -p {params.dir}
        
        /home/mattia/miniconda3/envs/samtools/bin/samtools sort -o {output.bam} -O bam {input.sam} -@ {threads}
        
        /home/mattia/miniconda3/envs/samtools/bin/samtools index {output.bam} -@ {threads}
        """


rule bam__dedup:
    input:
        markdup = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/{mysample}.bam"
    output:
        dedup = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam"),
        bai   = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam.bai"),
        metrics = temp("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam_metrics.txt")
    params:
        dir="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/",
        tmp = "/proj/tmp/tmp_MZ/tmp"
    resources:
        mem_mb=140000
    threads: 20
    #conda:
        #"picard"
    shell:
        """
        mkdir -p {params.dir}
        
        java -Xmx120g -jar /home/mattia/picard.jar MarkDuplicates I={input.markdup} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp}
        
        /home/mattia/miniconda3/envs/samtools/bin/samtools index {output.dedup} -@ {threads}
        
        """
#  Non canonical chromosomes, Y and X unmapped and mitochondrial are removed
rule bam__filter:
    """
    Filter out Non-canonical chromosomes, Y and X  and mitochondrial
    """
    input:
        dedup  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/{mysample}.bam"
    output:
        filter = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam",
        bai="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.bai"
    params:
        dir    = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter"
    resources:
        mem_mb=140000
    threads: 20
    #conda:
        #"samtools"
    shell:
        """
        mkdir -p {params.dir}

        /home/mattia/miniconda3/envs/samtools/bin/samtools view -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > {output.filter} -@ {threads}
        
        /home/mattia/miniconda3/envs/samtools/bin/samtools index {output.filter} -@ {threads}
        
        """
rule stat:
    """
    Filter out Non-canonical chromosomes, Y and X  and mitochondrial
    """
    input:
        filter = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam"
    output:
        flagstat = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.flagstat"
    resources:
        mem_mb=140000
    threads: 20
    #conda:
        #"samtools"
    shell:
        """        
        /home/mattia/miniconda3/envs/samtools/bin/samtools flagstat -@ {threads} {input.filter} > {output.flagstat} 
        
        """
#rule down_sample:
    #input:  
        #filter ="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam",
        #stat ="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam.flagstat"
    #output: 
        #downsample_bam="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam", 
        #downsample_bai="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam.bai"
    #resources: name = "downsample", time_min=50000, mem_mb=64000, cpus=20
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
        #else:
            #down_rate = 1

        #shell("/proj/tmp/tmp_MZ/anaconda3/envs/sambamba/bin/sambamba view -f bam -t 36 --subsampling-seed=3 -s {rate} {inbam} |  /proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} ".format(rate = down_rate, inbam = input[0], outbam = output[0]))

        #shell("/proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools index {outbam}".format(outbam = output[0]))

#rule downsample_stat:
    #input:
        #downsample_bam = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam"
    #output:
        #downsample_flagstat = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam.flagstat"
    #resources: name = "downsample_stat", time_min=50000, mem_mb=64000, cpus=20
    #shell:
        #"""        
        #/proj/tmp/tmp_MZ/anaconda3/envs/samtools/bin/samtools flagstat -@ 20 {input.downsample_bam} > {output.downsample_flagstat} 
        
        #"""

#rule bam__bigWig_downsample:
    #input: 
        #downsample_bam="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam",
        #downsample_bai="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}.bam.bai"
    #output:
        #downsample_bw="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample}_RPKM.bw"
    #params:
        #ref_genome_fa   =  config['ref_genome_fa'],
        # blacklist       =  config['blacklist'],
        #genome_size_bp  = config['genome_size_bp'],
        #mapping_qual_bw = 10
    #resources: name = "bamcoverage_downsample", time_min=50000, mem_mb=64000, cpus=10
    #shell:
        #"""
        #/proj/tmp/tmp_MZ/anaconda3/envs/deeptools/bin/bamCoverage -b {input.downsample_bam} --outFileName {output.downsample_bw} --normalizeUsing RPKM --binSize 10 --smoothLength 300 --numberOfProcessors 10 --effectiveGenomeSize {params.genome_size_bp}  --ignoreDuplicates  --skipNAs --exactScaling 
          
        #"""
#rule peak_calling_ATAC_downsample:
    #input:
        #downsample_treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mysample_atac}.bam",
    #output:
        #downsample_peaks_atac =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/SE/{mysample_atac}_peaks.narrowPeak",
        #downsample_bed   =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/SE/{mysample_atac}_summits.bed"
    #params:
        #peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/SE",
        #gsize      = config['genome_size_bp'],
        #qvalue     = config['peaks_qvalue'],
        #outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/SE"
    #resources: name = "macs2_ATAC_downsample", time_min=50000, mem_mb=64000, cpus=1
    #conda:
        #"macs2"
    #shell:
       #"""
        #mkdir -p {params.peaks_dir}
    
        #macs2 callpeak -t {input.downsample_treatment} --name {wildcards.mysample_atac} --outdir {params.outdir} --gsize {params.gsize} --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all  --qvalue {params.qvalue}  
    
        #"""

rule Gopeaks_P4_downsample:
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
    #conda:
        #"gopeaks"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        /home/mattia/miniconda3/envs/gopeaks/bin/gopeaks -b {input.downsample_treatment} -o {params.outdir}/{wildcards.pups} -p {params.qvalue} 
        
        """

rule Gopeaks_mature_downsample:
    input:
        downsample_treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/{mature}.bam",
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{ig2m}.bam",ig2m=IG2M, myrun=RUNID)
    output:
        downsample_peaks_mature =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/Gopeaks/{mature}_peaks.bed"
    params:
        peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/Gopeaks",
        qvalue     = config['peaks_qvalue'],
        outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/downsample/peaks/Gopeaks"
    resources:
        mem_mb=64000
    threads: 20
    #conda:
        #"/home/mattia/miniconda3/envs/gopeaks"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        /home/mattia/miniconda3/envs/gopeaks/bin/gopeaks -b {input.downsample_treatment} -o {params.outdir}/{wildcards.mature} -p {params.qvalue}
        
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
    #conda:
        #"deeptools"
    shell:
        """
        /home/mattia/miniconda3/envs/deeptools/bin/bamCoverage -b {input.bam} --outFileName {output.bw} --normalizeUsing RPKM --binSize 10 --smoothLength 300 --numberOfProcessors 10 --effectiveGenomeSize {params.genome_size_bp}  --ignoreDuplicates  --skipNAs --exactScaling  
          
        """

#rule bam__compare:
    #input: 
        #bam="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}.bam",
        #control=expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{control}.bam", control=CONTROL, myrun=RUNID)
    #output:
        #ratio="/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}_ratio_RPKM.bw"
    #params:
        #ref_genome_fa   =  config['ref_genome_fa'],
        # blacklist       =  config['blacklist'],
        #genome_size_bp  = config['genome_size_bp'],
        #mapping_qual_bw = 10
    #resources: name = "bamcompare", time_min=50000, mem_mb=64000, cpus=10
    #shell:
        #"""
        #/proj/tmp/tmp_MZ/anaconda3/envs/deeptools/bin/bamCompare -b1 {input.bam} -b2 {input.control} --outFileName {output.ratio} --normalizeUsing RPKM --binSize 10 --smoothLength 300 --numberOfProcessors 10 --scaleFactorsMethod None --effectiveGenomeSize {params.genome_size_bp}  --ignoreDuplicates  --skipNAs --exactScaling  
          
        #"""
rule peak_calling_ATAC:
    input:
        treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample_atac}.bam",
    output:
        peaks_atac =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/SE/{mysample_atac}_peaks.narrowPeak",
        bed   =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/SE/{mysample_atac}_summits.bed"
    params:
        peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/SE",
        gsize      = config['genome_size_bp'],
        qvalue     = config['peaks_qvalue'],
        outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/SE"
    resources: name = "macs2_ATAC", time_min=50000, mem_mb=64000, cpus=1
    #conda:
        #"macs2"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        /home/mattia/miniconda3/envs/macs2/bin/macs2 callpeak -t {input.treatment} --name {wildcards.mysample_atac} --outdir {params.outdir} --gsize {params.gsize} --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all  --qvalue {params.qvalue}  
    
        """

#rule peak_calling_narrow:
    #input:
         #treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample_narrow}.bam",
    #output:
        #peaks_narrow =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/{mysample_narrow}_peaks.narrowPeak",
        #bed   =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/{mysample_narrow}_summits.bed"
    #params:
        #peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/",
        #gsize      = config['genome_size_bp'],
        #qvalue     = config['peaks_qvalue'],
        #outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/"
    #resources: name = "macs2 Narrowpeak", time_min=50000, mem_mb=64000, cpus=20
    #shell:
        #"""
        #mkdir -p {params.peaks_dir}
    
        #/proj/tmp/tmp_MZ/anaconda3/envs/macs2/bin/macs2 callpeak -t {input.treatment} --name {wildcards.mysample_narrow} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --nomodel --qvalue {params.qvalue} --keep-dup all --call-summits
        
        #"""


#rule peak_calling_broad:
    #input:
        #treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample_broad}.bam",
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{input}.bam",myrun=RUNID,input=INPUT)
    #output:
        #peaks_broad =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/{mysample_broad}_peaks.broadPeak"
    #params:
        #peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/",
        #gsize      = config['genome_size_bp'],
        #qvalue     = config['peaks_qvalue'],
        #outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/"
    #resources: name = "macs2_broad", time_min=50000, mem_mb=64000, cpus=1
    #shell:
        #"""
        #mkdir -p {params.peaks_dir}
    
        #/proj/tmp/tmp_MZ/anaconda3/envs/macs2/bin/macs2 callpeak -t {input.treatment} -c {input.control} --name {wildcards.mysample_broad} --outdir {params.outdir} --gsize {params.gsize} -f BAMPE --broad --nomodel --qvalue {params.qvalue} --keep-dup all
        
        #"""
        
rule Gopeaks_mature:
    input:
        treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mature}.bam",
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{ig2m}.bam",ig2m=IG2M, myrun=RUNID)
    output:
        peaks_mature =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/Gopeaks/{mature}_peaks.bed"
    params:
        peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/Gopeaks",
        qvalue     = config['peaks_qvalue'],
        outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/Gopeaks"
    resources:
        mem_mb=64000
    threads: 20
    #conda:
        #"gopeaks"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        /home/mattia/miniconda3/envs/gopeaks/bin/gopeaks -b {input.treatment}  -o {params.outdir}/{wildcards.mature} -p {params.qvalue} --broad
        
        """
        
rule Gopeaks_P4:
    input:
        treatment = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{pups}.bam",
        #control= expand("/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{igp4}.bam",igp4=IGP4, myrun=RUNID)
    output:
        peaks_pups =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/Gopeaks/{pups}_peaks.bed"
    params:
        peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/Gopeaks",
        qvalue     = config['peaks_qvalue'],
        outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/Gopeaks"
    resources:
        mem_mb=64000
    threads: 20
    #conda:
        #"/home/mattia/miniconda3/envs/gopeaks"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        /home/mattia/miniconda3/envs/gopeaks/bin/gopeaks -b {input.treatment} -o {params.outdir}/{wildcards.pups} -p {params.qvalue} --broad
        
        """
        
#rule merge:
    #input:
        #peaks = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/{mysample}_peaks.narrowPeak",
    #output:
        #cat =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/{mysample_broad}_peaks.broadPeak",
        #bed   =  "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/{mysample_broad}_summits.bed"
    #params:
        #peaks_dir  = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/",
        #gsize      = config['genome_size_bp'],
        #qvalue     = config['peaks_qvalue'],
        #outdir     = "/proj/tmp/tmp_MZ/{myrun}/fastq/trimmed/trimmomatic/mapped/bowtie2/sorted/samtools/dedup/filter/peaks/"
    #shell:
        #"""
        #mkdir -p {params.peaks_dir}
        
        #cat {input.peaks}  > {output.merged}
        
        #/proj/tmp/tmp_MZ/anaconda3/envs/bedtools/envs/bedtools sort  {output.merged} > {output.merged}
        
        #/proj/tmp/tmp_MZ/anaconda3/envs/bedtools/envs/bedtools merge  {output.merged} > {output.merged}
        
        #"""
