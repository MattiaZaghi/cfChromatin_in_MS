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
TAG= config["tag"]
CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
CHIPS_SE = [sample for sample in MARK_SAMPLES if CHIP in sample]
TAGS=[sample for sample in MARK_SAMPLES if TAG in sample]
RUNID = config["RUN_ID"]
## list BAM files
ALL_SAMPLES = TAGS
BEDGRAPH=expand("{myrun}/bedgraph/bedtools/{sample}.bedgraph", sample=ALL_SAMPLES, myrun=RUNID)
ALL_BIGWIG= expand("{myrun}/coverage/bedgraphtobigwig/{sample}_CPM.bw", sample = ALL_SAMPLES,myrun=RUNID)
#MACS2 = expand("{myrun}/peaks/macs2/{sample}.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
#MACS2_BROAD = expand("{myrun}/peaks/macs2/{sample}.broadPeak", sample = ALL_SAMPLES,myrun=RUNID)



TARGETS = []
TARGETS.extend(BEDGRAPH)
TARGETS.extend(ALL_BIGWIG)
#TARGETS.extend(MACS2)
#TARGETS.extend(MACS2_BROAD)






ruleorder:  unzip_tag > tag_to_bedgraph > coverage > macs2 > macs2_broad 




rule all:
    input: TARGETS

rule unzip_tag:
    input:
        tag_gz=lambda wildcards: FILES[wildcards.sample.split('_')[0]][wildcards.sample.split('_')[1]][wildcards.sample.split('_')[2]][0]
    output:
        tag="{myrun}/cat/{sample}.tagAlign"
    threads: config['THREADS']
    shell:
        """
        gunzip -c {input.tag_gz} > {output.tag}

        """


rule tag_to_bedgraph:
    """
    Align reads with bowtie2
    """
    input:
        tag="{myrun}/cat/{sample}.tagAlign"
    output:
        bedgraph = "{myrun}/bedgraph/bedtools/{sample}.bedgraph",
        sorted = "{myrun}/bedgraph/bedtools/sorted/{sample}.bedgraph"
    params:
        dir      = "{myrun}/mapped/bedgraph",
        chrom    = config['chrom_size']
    threads: config['THREADS']
    resources:
        mem_mb=64000
    conda:
        "/home/mattia/miniconda3/envs/bedtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        bedtools genomecov -i {input.tag} -bg  -g {params.chrom} > {output.bedgraph}

        bedtools sort -i {output.bedgraph} > {output.sorted}

        """

rule coverage:
    input: 
        sorted = "{myrun}/bedgraph/bedtools/sorted/{sample}.bedgraph"
    output:
        bw="{myrun}/coverage/bedgraphtobigwig/{sample}_CPM.bw"
    params:
        chrom    = config['chrom_size']
    resources:
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/coverage/bedgraphtobigwig/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/ucsc.yml"
    shell:
        """
        bedGraphToBigWig {input.sorted}  {params.chrom} {output.bw}

        """

rule macs2:
    input:
        sorted = "{myrun}/bedgraph/bedtools/sorted/{sample}.bedgraph"
    output:
        peaks_narrow =  "{myrun}/peaks/macs2/{sample}.narrowPeak"
    params:
        peaks_dir  = "{myrun}/peaks/macs2/"
    resources: mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/macs2.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        macs2 bdgpeakcall -i {input.sorted} -c 1.0 -l 100 -g 500 -o {output.peaks_narrow}
        
        """

rule macs2_broad:
    input:
        sorted = "{myrun}/bedgraph/bedtools/sorted/{sample}.bedgraph"
    output:
        peaks_broad =  "{myrun}/peaks/macs2/{sample}.broadPeak"
    params:
        peaks_dir  = "{myrun}/peaks/macs2/"
    resources: mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/macs2.yml"
    shell:
        """
        mkdir -p {params.peaks_dir}
    
        macs2 bdgbroadcall -i {input.sorted} -o {output.peaks_broad} --cutoff-peak 1.0 --minlen 500 --maxgap 1000 --cutofflink 1.5 --maxgaplink 500 
        
        """
        
        