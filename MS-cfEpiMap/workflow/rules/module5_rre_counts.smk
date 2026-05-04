"""
Module 5 — RRE Count Matrices
==============================
Counts raw fragment midpoints at each RRE subset per sample.
Produces integer count matrices (regions × samples) fed directly into DESeq2.
"""

# Map subset name → BED file path.
# Per-cell-type CNS subsets (cns_oligodendrocyte, cns_neuron, …) are added
# dynamically from config so the pipeline automatically runs DESeq2 on each.
RRE_BED = {
    "full":          config["rre"]["universe_bed"],
    "cns":           config["rre"]["cns_bed"],
    "immune":        config["rre"]["immune_bed"],
    "gwas_proximal": config["rre"]["gwas_proximal_bed"],
    **{
        f"cns_{ct}": config["rre"]["cns_celltype_beds"][ct]
        for ct in config["rre"].get("cns_cell_types", [])
    },
}

rule count_rre_midpoints_per_sample:
    """Count midpoints overlapping one RRE subset for one sample."""
    input:
        frags   = lambda wc: f"{config['data']['frags_dir']}/{wc.sample}.bed",
        regions = lambda wc: RRE_BED[wc.subset],
    output:
        temp("results/counts/rre_per_sample/{subset}/{sample}.counts"),
    params:
        script = workflow.basedir + "/scripts/count_midpoints.py",
    threads: 1
    conda:  "../../envs/python_analysis.yaml"
    shell:
        """
        python {params.script} \
            --frags   {input.frags} \
            --regions {input.regions} \
            --output  {output}
        """

rule build_rre_count_matrix:
    """Assemble per-sample counts into a regions × samples count matrix."""
    input:
        counts  = lambda wc: expand(
            "results/counts/rre_per_sample/{subset}/{sample}.counts",
            subset=wc.subset,
            sample=SAMPLES,
        ),
        regions = lambda wc: RRE_BED[wc.subset],
    output:
        "results/counts/rre_{subset}_counts.tsv",
    params:
        samples = SAMPLES,
        subset  = "{subset}",
    conda:  "../../envs/python_analysis.yaml"
    script: "../../scripts/module5_rre_counts.py"
