"""
External H3K27ac cfChIP-seq enrichment score module.

Computes per-sample enrichment scores for published fragment BED files
(Sadeh et al., GSE accession) and cross-compares with the internal MS cohort.

Rules:
  compute_external_enrichment_score   — per GSM sample
  aggregate_external_enrichment_scores — merge all per-sample TSVs
  enrichment_score_cross_comparison    — internal vs external combined plot
"""

import glob as _glob
import os as _os

_EXT_BED_DIR = config.get('external_bed_dir', '/date/gcb/gcb_MZ/Analysis/BED/H3K27ac')
_ENRICH_ON   = config.get('chromhmm_ontarget',  'ref_files/chromhmm_ontarget_H3K27ac.bed')
_ENRICH_OFF  = config.get('chromhmm_offtarget', 'ref_files/chromhmm_offtarget_H3K27ac.bed')
_OUT         = config['outputFolder']

# Auto-discover all GSM*.bed files in the external directory
_EXT_BED_FILES = sorted(_glob.glob(_os.path.join(_EXT_BED_DIR, 'GSM*.bed')))
_EXT_SAMPLES   = [_os.path.splitext(_os.path.basename(f))[0] for f in _EXT_BED_FILES]

if not _EXT_SAMPLES:
    print(f"Warning: no GSM*.bed files found in {_EXT_BED_DIR}; "
          "external enrichment targets will be skipped")


rule compute_external_enrichment_score:
    """
    Compute H3K27ac enrichment score for one external (GSM) sample.
    Reuses compute_enrichment_score.py directly on the fragment BED —
    no prior pipeline steps are needed for these published BED files.
    """
    conda: "envs/qc.yaml"
    input:
        frag_bed  = lambda wc: _os.path.join(_EXT_BED_DIR, wc.ext_sample + '.bed'),
        ontarget  = _ENRICH_ON,
        offtarget = _ENRICH_OFF,
        script    = "auxiliar_programs/compute_enrichment_score.py"
    output:
        tsv = _OUT + "/external_enrichment/{ext_sample}_enrichment.tsv"
    params:
        sample = "{ext_sample}"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        python {input.script} \
            --frag_bed  {input.frag_bed} \
            --ontarget  {input.ontarget} \
            --offtarget {input.offtarget} \
            --sample    {params.sample} \
            --out       {output.tsv}
        """


rule aggregate_external_enrichment_scores:
    """
    Merge all per-GSM-sample enrichment TSVs into a single table.
    Output: reports/enrichment/external_enrichment_scores.tsv
    """
    input:
        tsvs = expand(
            _OUT + "/external_enrichment/{ext_sample}_enrichment.tsv",
            ext_sample = _EXT_SAMPLES
        )
    output:
        aggregate = _OUT + "/reports/enrichment/external_enrichment_scores.tsv"
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


rule enrichment_score_cross_comparison:
    """
    Combined violin + boxplot comparing internal MS cohort enrichment scores
    against the external published Sadeh et al. cohort.
    Groups ordered by median enrichment score; internal = blues,
    external healthy = green, external cancer = reds/oranges.
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        internal = _OUT + "/reports/enrichment/all_enrichment_scores.tsv",
        external = _OUT + "/reports/enrichment/external_enrichment_scores.tsv",
        script   = "auxiliar_programs/enrichment_score_cross_comparison.R"
    output:
        pdf = _OUT + "/reports/enrichment/enrichment_cross_comparison.pdf"
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        Rscript {input.script} \
            {input.internal} \
            {input.external} \
            {output.pdf} \
            2>&1 || true
        touch {output.pdf}
        """
