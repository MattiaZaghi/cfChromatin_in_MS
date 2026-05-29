# Base directory containing the BED files
BED_DIR = "/date/gcb/gcb_MZ/Analysis/BED/H3K27ac"

# Discover basenames of all *.bed files in BED_DIR
SAMPLES = sorted(glob_wildcards(f"{BED_DIR}" + "/{sample}.bed").sample)  # e.g. sample = filename without .bed [web:23][web:32]

# Define final targets
rule all:
    input:
        expand(f"{BED_DIR}" + "/{{sample}}.future_yield.txt", sample=SAMPLES)  # create one output per input BED [web:23][web:24]


# Run preseq on each BED (use sorted if produced)
rule preseq_lc_extrap:
    input:
        lambda wildcards: (
            f"{BED_DIR}/{wildcards.sample}.bed"
            if os.path.exists(f"{BED_DIR}/{wildcards.sample}.sorted.bed")
            else f"{BED_DIR}/{wildcards.sample}.bed"
        )  # choose sorted if present [web:24]
    output:
        f"{BED_DIR}" + "/{sample}.future_yield.txt"  # unique per sample [web:27]
    threads: 1
    conda:
        "/home/mattia/miniconda3/envs/preseq.yml"
    resources:
        mem_mb=2000
    shell:
        # For BED input, do not use -B; for BAM, switch input and add -B [web:27]
        "preseq lc_extrap -o {output} {input} "
