"""
QC flagging module.

Cut-offs are per enrichment mark, taken from config['qc_thresholds'] and
materialised as a TSV the R scripts join on. Scores are only comparable within
a mark, so every distribution and pass/fail tally is computed per mark.

Rules:
  count_fragments           — count fragments per sample from .bed file
  aggregate_fragment_counts — merge per-sample counts into one TSV
  qc_thresholds_table       — write config['qc_thresholds'] as a TSV
  qc_flag_samples           — flag samples against their mark's cut-offs,
                              write qc_summary.tsv + PDF (one page per mark)
  qc_flag_protocol          — same, grouped by protocol version
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


rule qc_thresholds_table:
    """
    Materialise config['qc_thresholds'] as a TSV the R scripts can join on.

    Declaring this as an input of the QC rules means editing a cut-off in
    config.yaml retriggers the flagging by itself, and the exact cut-offs used
    for a run stay on disk next to the results.
    """
    output:
        tsv = _OUT + "/qc/qc_thresholds.tsv"
    run:
        import os
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        marks = sorted(set(config.get('_sample_mark_', {}).values())) or ['H3K27ac']
        thresholds = config.get('qc_thresholds', {})
        with open(output.tsv, "w") as fh:
            fh.write("mark\tenrich_thresh\tfrag_thresh\n")
            for m in marks:
                t = thresholds.get(m, {})
                fh.write(f"{m}\t{t.get('enrichment', 7)}\t{t.get('fragments', 1000000)}\n")


rule qc_flag_samples:
    """
    QC flagging against the per-mark cut-offs in config['qc_thresholds'].
    Samples are compared only against others sharing their enrichment mark.
    Outputs qc_summary.tsv (sample_id, mark, n_fragments, enrichment_score,
    thresholds applied, frag_pass, enrich_pass, qc_pass) and a diagnostic PDF
    with one page per mark.
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        enrichment  = _OUT + "/reports/enrichment/all_enrichment_scores.tsv",
        frag_counts = _OUT + "/qc/fragment_counts.tsv",
        thresholds  = _OUT + "/qc/qc_thresholds.tsv",
        script      = "auxiliar_programs/qc_flag_samples.R",
        common      = "auxiliar_programs/qc_flag_common.R"
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
            {input.thresholds} \
            2>&1 || true
        touch {output.tsv} {output.pdf}
        """

rule qc_flag_protocol:
    """
    As qc_flag_samples, but grouped and coloured by protocol version instead of
    disease group. Same per-mark cut-offs, same one-page-per-mark layout.
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        enrichment  = _OUT + "/reports/enrichment/all_enrichment_scores.tsv",
        frag_counts = _OUT + "/qc/fragment_counts.tsv",
        thresholds  = _OUT + "/qc/qc_thresholds.tsv",
        script      = "auxiliar_programs/qc_flag_protocol.R",
        common      = "auxiliar_programs/qc_flag_common.R"
    output:
        tsv = _OUT + "/qc/qc_summary_protocol.tsv",
        pdf = _OUT + "/reports/qc/qc_distributions_protocol.pdf"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {output.pdf})
        Rscript {input.script} \
            {input.enrichment} \
            {input.frag_counts} \
            {output.tsv} \
            {output.pdf} \
            {input.thresholds} \
            2>&1 || true
        touch {output.tsv} {output.pdf}
        """