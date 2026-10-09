"""
cfChIP-seq enrichment score module.

Computes per-sample enrichment = (on-target density) / (off-target density)
where on/off-target regions come from the 18-state Epimap chromHMM
segmentation (build_chromhmm_reference.py, run once before the pipeline).

The reference pair is chosen per sample from the samplesheet's
'enrichment_mark' column, so H3K4me3 libraries are scored against the H3K4me3
region set and H3K27ac libraries against the H3K27ac one.

Rules:
  compute_enrichment_score   — per sample, uses fragment midpoints
  aggregate_enrichment_scores — merge all per-sample TSVs into one table
  enrichment_score_plot       — per-sample ECDF plot with knee-point
"""

_OUT = config['outputFolder']


rule compute_enrichment_score:
    """
    Compute the enrichment score for one sample against its own mark's regions.
    Score = (on_frags / on_bp) / (off_frags / off_bp)
    Uses fragment midpoints (bedtools intersect) against chromHMM BEDs.
    """
    conda: "envs/qc.yaml"
    input:
        frag_bed  = _OUT + "/frags/{sample}.bed",
        ontarget  = lambda wc: enrichment_refs(sample_mark(wc.sample))[0],
        offtarget = lambda wc: enrichment_refs(sample_mark(wc.sample))[1],
        script    = "auxiliar_programs/compute_enrichment_score.py"
    output:
        tsv = _OUT + "/enrichment/{sample}_enrichment.tsv"
    params:
        sample = "{sample}",
        mark   = lambda wc: sample_mark(wc.sample)
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        python {input.script} \
            --frag_bed  {input.frag_bed} \
            --ontarget  {input.ontarget} \
            --offtarget {input.offtarget} \
            --sample    {params.sample} \
            --mark      {params.mark} \
            --out       {output.tsv}
        """


rule aggregate_enrichment_scores:
    """
    Merge all per-sample enrichment TSVs into one table.
    Output: reports/enrichment/all_enrichment_scores.tsv
    """
    input:
        tsvs = expand(
            _OUT + "/enrichment/{sample}_enrichment.tsv",
            sample = config.get('_samples_', [])
        )
    output:
        aggregate = _OUT + "/reports/enrichment/all_enrichment_scores.tsv"
    run:
        import os
        os.makedirs(os.path.dirname(output.aggregate), exist_ok=True)
        with open(output.aggregate, "w") as out_fh:
            header_written = False
            for tsv in input.tsvs:
                with open(tsv) as fh:
                    for i, line in enumerate(fh):
                        if i == 0:
                            if not header_written:
                                out_fh.write(line)
                                header_written = True
                        else:
                            out_fh.write(line)


rule enrichment_score_comparison_plots:
    """
    Cohort-level comparison boxplots of enrichment scores, one page per mark:
      - By protocol version (V1/V2/V3/V4/1D — inferred from sample IDs)
      - By sample type (Plasma / CSF — inferred from sample IDs)
      - By disease group (Ctrl / MS-New / MS-Rit-Stable / …)
    Kruskal-Wallis p-value annotated on each panel, computed within a mark.
    Output: reports/enrichment/enrichment_comparison_plots.pdf
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        aggregate = _OUT + "/reports/enrichment/all_enrichment_scores.tsv",
        script    = "auxiliar_programs/enrichment_score_comparison_plots.R"
    output:
        pdf = _OUT + "/reports/enrichment/enrichment_comparison_plots.pdf"
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        Rscript {input.script} \
            {input.aggregate} \
            {output.pdf} \
            2>&1 || true
        touch {output.pdf}
        """


rule enrichment_score_plot:
    """
    Per-sample enrichment score QC plot.
    Shows the ECDF of all samples sharing this sample's mark, with:
      - current sample highlighted (red dot + label)
      - quantile dotted lines Q10/Q25/Q50/Q75/Q90
      - knee-point of the distribution marked (green line)
    X-axis: enrichment score value (optimized score cutoff)
    Y-axis: cumulative fraction of samples
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        aggregate = _OUT + "/reports/enrichment/all_enrichment_scores.tsv",
        script    = "auxiliar_programs/enrichment_score_plot.R"
    output:
        pdf = _OUT + "/reports/enrichment/{sample}_enrichment_score.pdf"
    params:
        sample = "{sample}"
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        Rscript {input.script} \
            {input.aggregate} \
            {params.sample} \
            {output.pdf} \
            2>&1 || true
        touch {output.pdf}
        """
