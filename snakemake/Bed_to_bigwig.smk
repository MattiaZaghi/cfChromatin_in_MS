import os
from glob import glob

configfile: "/home/mattia/cfChromatin_in_MS/snakemake/config_ChIP.yaml"

CHROM_SIZES = config["chrom_sizes"]
OUTDIR = config.get("outdir", "results")
INPUT_GLOB = config.get("input_beds_glob", "/proj/user/mattia/healthy_baca/*.bed")
NORMALIZE_CPM = bool(config.get("normalize_cpm", True))
THREADS = int(config.get("threads", 4))
BLACKLIST = config.get("blacklist_bed", "")
STRAND = config.get("strand", "none")

# Discover inputs
BEDS = sorted(glob(INPUT_GLOB))
SAMPLES = [os.path.splitext(os.path.basename(x))[0] for x in BEDS]

# Helper for input (OK), not for output
def outpath(sample, ext):
    return os.path.join(OUTDIR, sample + ext)

# Sanity: if SAMPLES is empty, warn early
if len(SAMPLES) == 0:
    print(f"[WARN] No BED files found for pattern: {INPUT_GLOB}")
    print("       Update 'input_beds_glob' in config to point to your BED files.")

rule all:
    input:
        # Ensure mkdirs runs by depending on OUTDIR as a directory
        expand(os.path.join(OUTDIR, "{sample}.raw.bw"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}.cpm.bw"), sample=SAMPLES) if NORMALIZE_CPM else []


rule sort_bed:
    input:
        # Dynamic input function allowed here
        bed=lambda wc: next(x for x in BEDS if os.path.splitext(os.path.basename(x))[0] == wc.sample),
    output:
        bed=os.path.join(OUTDIR, "{sample}.sorted.bed")
    threads: THREADS
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    shell:
        r"""
        set -euo pipefail
        LC_ALL=C sort -k1,1 -k2,2n "{input.bed}" > "{output.bed}"
        """

rule filter_blacklist:
    input:
        bed=os.path.join(OUTDIR, "{sample}.sorted.bed")
    output:
        bed=os.path.join(OUTDIR, "{sample}.sorted.filt.bed")
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    run:
        import shutil
        if not BLACKLIST:
            shutil.copy(input.bed, output.bed)
        else:
            shell(r"""
                bedtools intersect -v -a "{input.bed}" -b "{BLACKLIST}" > "{output.bed}"
            """)

rule bed_to_bedgraph_raw:
    input:
        bed=os.path.join(OUTDIR, "{sample}.sorted.filt.bed"),
        chrom=CHROM_SIZES
    output:
        bg=os.path.join(OUTDIR, "{sample}.raw.bedGraph")
    params:
        strand=STRAND
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    shell:
        r"""
        set -euo pipefail

        if [ "{params.strand}" = "none" ]; then
            bedtools genomecov -i "{input.bed}" -g "{input.chrom}" -bg > "{output.bg}"
        else
            bedtools genomecov -i "{input.bed}" -g "{input.chrom}" -bg -strand "{params.strand}" > "{output.bg}"
        fi

        # Ensure sorting for bedGraphToBigWig
        LC_ALL=C sort -k1,1 -k2,2n -o "{output.bg}" "{output.bg}"
        """

rule bedgraph_to_bigwig_raw:
    input:
        bg=os.path.join(OUTDIR, "{sample}.raw.bedGraph"),
        chrom=CHROM_SIZES
    output:
        bw=os.path.join(OUTDIR, "{sample}.raw.bw")
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    shell:
        r"""
        set -euo pipefail
        bedGraphToBigWig "{input.bg}" "{input.chrom}" "{output.bw}"
        """

rule count_reads:
    input:
        bed=os.path.join(OUTDIR, "{sample}.sorted.filt.bed")
    output:
        count_file=os.path.join(OUTDIR, "{sample}.count.txt")
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    shell:
        r"""
        set -euo pipefail
        wc -l < "{input.bed}" | awk '{{print $1}}' > "{output.count_file}"
        """
rule cpm_scale_bedgraph:
    input:
        bg=os.path.join(OUTDIR, "{sample}.raw.bedGraph"),
        count_file=os.path.join(OUTDIR, "{sample}.count.txt")
    output:
        bg=os.path.join(OUTDIR, "{sample}.cpm.bedGraph")
    params:
        scale=lambda wildcards, input: 1e6 / int(open(input.count_file).read().strip() or 1)
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    shell:
        r"""
        awk -v OFS='\t' -v s={params.scale} '{{{{print $1,$2,$3,$4*s}}}}' "{input.bg}" > "{output.bg}"
        LC_ALL=C sort -k1,1 -k2,2n -o "{output.bg}" "{output.bg}"
        """

rule bedgraph_to_bigwig_cpm:
    input:
        bg=os.path.join(OUTDIR, "{sample}.cpm.bedGraph"),
        chrom=CHROM_SIZES
    output:
        bw=os.path.join(OUTDIR, "{sample}.cpm.bw")
    conda:
        "/home/mattia/miniconda3/envs/bed_to_bigwig.yaml"
    shell:
        r"""
        set -euo pipefail
        bedGraphToBigWig "{input.bg}" "{input.chrom}" "{output.bw}"
        """

# Optional: if CPM disabled, all targets are raw only
if not NORMALIZE_CPM:
    ruleorder: bedgraph_to_bigwig_raw > bedgraph_to_bigwig_cpm
