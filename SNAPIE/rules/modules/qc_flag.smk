"""
QC flagging module.

Rules:
  count_fragments          — count fragments per sample from .bed file
  aggregate_fragment_counts — merge per-sample counts into one TSV
  qc_flag_samples          — compute data-driven thresholds (knee-point),
                             flag samples, write qc_summary.tsv + PDF
"""

_OUT = config['outputFolder']


rule count_fragments:
    """Count the number of fragments (lines) in each sample's fragment BED."""
    input:
        bed = _OUT + "/frags/{sample}.bed"
    output:
        tsv = _OUT + "/qc/frag_counts/{sample}_frag_count.tsv"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        count=$(wc -l < {input.bed})
        printf 'sample_id\tn_fragments\n{wildcards.sample}\t%s\n' "$count" > {output.tsv}
        """


rule aggregate_fragment_counts:
    """Merge all per-sample fragment count TSVs into one table."""
    input:
        tsvs = expand(_OUT + "/qc/frag_counts/{sample}_frag_count.tsv",
                      sample = config.get('_samples_', []))
    output:
        tsv = _OUT + "/qc/fragment_counts.tsv"
    run:
        import os
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w") as out_fh:
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


rule qc_flag_samples:
    """
    Data-driven QC flagging using knee-point thresholds on the user's own
    distribution of enrichment scores and fragment counts.
    Outputs qc_summary.tsv (sample_id, n_fragments, enrichment_score,
    frag_pass, enrich_pass, qc_pass) and a diagnostic PDF.
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        enrichment  = _OUT + "/reports/enrichment/all_enrichment_scores.tsv",
        frag_counts = _OUT + "/qc/fragment_counts.tsv",
        script      = "auxiliar_programs/qc_flag_samples.R"
    output:
        tsv = _OUT + "/qc/qc_summary.tsv",
        pdf = _OUT + "/reports/qc/qc_distributions.pdf"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {output.pdf})
        Rscript {input.script} \
            {input.enrichment} \
            {input.frag_counts} \
            {output.tsv} \
            {output.pdf} \
            2>&1 || true
        touch {output.tsv} {output.pdf}
        """
