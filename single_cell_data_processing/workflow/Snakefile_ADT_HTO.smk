import yaml, os,json
cfg = config                              # Snakemake’s automatic dict

SAMPLES = list(cfg["samples"].keys())
KINDS   = ["adt", "hto"]

rule all:
    input:
        expand("{outdir}/{sample}/{kind}/alevin/quants_mat.gz",
               outdir = cfg["outdir"],
               sample = SAMPLES,
               kind   = KINDS)

###############################################################################
# Quantify FASTQs with salmon alevin 1.4.0 (no index building)
###############################################################################
rule alevin_quant:
    input:
        idx = lambda wc: cfg["indices"][wc.kind],                # <-- fixed
        r1  = lambda wc: cfg["samples"][wc.sample][f"{wc.kind}_r1"],
        r2  = lambda wc: cfg["samples"][wc.sample][f"{wc.kind}_r2"]
    output:
        quant = temp("{outdir}/{sample}/{kind}/alevin/quants_mat.gz")
    params:
        fs   = cfg["geometry"]["featureStart"],
        fl   = cfg["geometry"]["featureLength"],
        ec   = cfg["params"]["expectCells"],
        out  = lambda wc: f"{cfg['outdir']}/{wc.sample}/{wc.kind}"
    threads: cfg["params"]["threads"]
    conda: "/home/mattia/miniconda3/envs/salmon.yml"
    shell:
        """
        salmon alevin -l ISR -i {input.idx} \
              -1 {input.r1} -2 {input.r2} \
              -p {threads}  --chromium \
              --citeseq \
              --featureStart {params.fs} --featureLength {params.fl} \
              --expectCells {params.ec} --naiveEqclass \
              -o {params.out}

        touch {output.quant}
        """
