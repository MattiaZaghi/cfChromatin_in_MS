"""
Module 4 — Cell-Type Deconvolution
====================================
Estimates per-sample cell-type contributions via Poisson z-scores at
cell-type-specific H3K27ac signature regions.  Derives CDI, NII, IAI.

Two execution paths, selected by the config key `deconv_chunk_size`:

  chunk_size = 0 (DEFAULT)
      One `compute_deconvolution` job does everything, exactly as before.
      Nothing about existing runs changes.

  chunk_size = N > 0
      The rule is split, because the two halves of module4_deconvolution.py
      have very different memory profiles:

        cross  — cohort Poisson/NB z vs the Ctrl distribution.  Counts are
                 per-cell-type TOTALS (a handful of ints per sample), so
                 memory is flat in sample count.  Runs ONCE over ALL samples;
                 it has to, because the scaling factor is a cohort median and
                 the z-scores are relative to the Ctrl mean.  This is what
                 keeps signature_scores / composite_indices / the heatmap and
                 violin cohort-global rather than chunk-local.

        within — per-sample NB-GLM vs Sadeh local background.  Counts are one
                 int per region per sample and `_ws_sample_counts` holds every
                 sample at once, so memory is LINEAR in sample count.  This is
                 the half that OOMs (61.5 GiB at 990 samples against a 64000
                 mem_mb request) and the only half that needs chunking.

      The within half is numerically independent per sample -- the background
      model, alpha_s and the BH-FDR are all fit inside one sample, and it never
      touches the cohort-median scaling factor -- so scoring a chunk and
      concatenating gives byte-identical rows to one big job.

      NOT reproduced in the chunked path: drivers_dir / drivers_manifest /
      drivers_enrichment_pdf aggregate ACROSS samples, so the chunk rule simply
      does not declare them (module4_deconvolution.py skips them silently when
      they are absent).  The three within-sample figures are emitted per chunk
      and are chunk-local.  Use chunk_size = 0 when you need those artefacts.
"""

# Read as a FLAT top-level key, not config["deconvolution"]["chunk_size"].
# A `deconvolution:` block in a --configfile override would replace the one from
# the base config unless the merge is deep, taking off_target_bed with it; a
# flat key cannot do that. The nested form is still honoured as a fallback for
# anyone who sets it directly in the base config.
_DECONV_CHUNK = int(
    config.get("deconv_chunk_size",
               config.get("deconvolution", {}).get("chunk_size", 0)) or 0
)
_DECONV_CROSS_MEM  = int(config.get("deconv_cross_mem_mb", 32000))
_DECONV_WITHIN_MEM = int(config.get("deconv_within_mem_mb", 64000))


if _DECONV_CHUNK <= 0:

    rule compute_deconvolution:
        """
        For each cell type, count midpoints at cell-type-specific signature
        regions (present in that type, absent in all others).  Compute
        Poisson z-scores vs Ctrl distribution.  Derive composite indices.
        Run within-sample NB-GLM enrichment vs Sadeh-style local background
        built from a curated off-target BED.
        """
        input:
            frags_list     = [SAMPLE_FRAG_BED[s] for s in SAMPLES],
            rre_universe   = config["rre"]["universe_bed"],
            anchors        = config["normalization"]["constitutive_anchors_bed"],
            off_target_bed = config["deconvolution"]["off_target_bed"],   # NEW
            sf             = "results/normalization/constitutive_scaling_factors.tsv",
            anchor_matrix  = "results/normalization/anchor_counts.tsv",
            meta           = config["sample_metadata"],
        output:
            scores                     = "results/deconvolution/signature_scores.tsv",
            composite                  = "results/deconvolution/composite_indices.tsv",
            heatmap                    = "results/deconvolution/deconvolution_heatmap.pdf",
            violin                     = "results/deconvolution/deconvolution_violin.pdf",
            within_sample_scores       = "results/deconvolution/scores_within_sample.tsv",
            within_sample_qc           = "results/deconvolution/qc_within_sample.tsv",
            within_sample_log2fe       = "results/deconvolution/fig_within_sample_log2FE.pdf",
            within_sample_contribution = "results/deconvolution/fig_within_sample_contribution.pdf",
            within_sample_qc_png       = "results/deconvolution/qc_within_sample.png",
            within_sample_windows_bed  = "results/deconvolution/background_signal_windows.bed.gz",
            drivers_dir              = directory("results/deconvolution/drivers"),
            drivers_manifest         = "results/deconvolution/drivers_manifest.tsv",
            drivers_enrichment_pdf   = "results/deconvolution/drivers_enrichment.pdf",
        params:
            samples        = SAMPLES,
            sample_bed_map = SAMPLE_FRAG_BED,
            genome_fasta   = config["genome"].get("fasta", ""),
            gene_annotation_bed = config["deconvolution"].get("gene_annotation_bed", ""),
            stage          = "all",
        conda:   "../../envs/python_analysis.yaml"
        threads: config["threads"]
        resources:
            mem_mb = 64000,
        script:  "../../scripts/module4_deconvolution.py"


else:

    # Deterministic contiguous chunks; the split only has to be a partition of
    # SAMPLES, since every within-sample result depends on that sample alone.
    _DECONV_CHUNKS = {
        f"{_k:03d}": SAMPLES[_i:_i + _DECONV_CHUNK]
        for _k, _i in enumerate(range(0, len(SAMPLES), _DECONV_CHUNK))
    }
    _DECONV_CHUNK_IDS = sorted(_DECONV_CHUNKS)

    rule deconvolve_cross_sample:
        """Cohort-level half: every sample in one job, flat memory."""
        input:
            frags_list     = [SAMPLE_FRAG_BED[s] for s in SAMPLES],
            rre_universe   = config["rre"]["universe_bed"],
            anchors        = config["normalization"]["constitutive_anchors_bed"],
            off_target_bed = config["deconvolution"]["off_target_bed"],
            sf             = "results/normalization/constitutive_scaling_factors.tsv",
            anchor_matrix  = "results/normalization/anchor_counts.tsv",
            meta           = config["sample_metadata"],
        output:
            scores    = "results/deconvolution/signature_scores.tsv",
            composite = "results/deconvolution/composite_indices.tsv",
            heatmap   = "results/deconvolution/deconvolution_heatmap.pdf",
            violin    = "results/deconvolution/deconvolution_violin.pdf",
        params:
            samples        = SAMPLES,
            sample_bed_map = SAMPLE_FRAG_BED,
            genome_fasta   = config["genome"].get("fasta", ""),
            gene_annotation_bed = config["deconvolution"].get("gene_annotation_bed", ""),
            stage          = "cross",
        conda:   "../../envs/python_analysis.yaml"
        threads: config["threads"]
        resources:
            mem_mb = _DECONV_CROSS_MEM,
        script:  "../../scripts/module4_deconvolution.py"

    rule deconvolve_within_chunk:
        """Per-sample half, one job per chunk of samples."""
        input:
            frags_list     = lambda w: [SAMPLE_FRAG_BED[s] for s in _DECONV_CHUNKS[w.chunk]],
            rre_universe   = config["rre"]["universe_bed"],
            anchors        = config["normalization"]["constitutive_anchors_bed"],
            off_target_bed = config["deconvolution"]["off_target_bed"],
            sf             = "results/normalization/constitutive_scaling_factors.tsv",
            anchor_matrix  = "results/normalization/anchor_counts.tsv",
            meta           = config["sample_metadata"],
        output:
            within_sample_scores       = temp("results/deconvolution/chunks/{chunk}/scores_within_sample.tsv"),
            within_sample_qc           = temp("results/deconvolution/chunks/{chunk}/qc_within_sample.tsv"),
            within_sample_windows_bed  = temp("results/deconvolution/chunks/{chunk}/background_signal_windows.bed.gz"),
            # Chunk-local figures: emitted so the script has somewhere to write,
            # not used downstream. Cohort versions need chunk_size = 0.
            within_sample_log2fe       = temp("results/deconvolution/chunks/{chunk}/fig_within_sample_log2FE.pdf"),
            within_sample_contribution = temp("results/deconvolution/chunks/{chunk}/fig_within_sample_contribution.pdf"),
            within_sample_qc_png       = temp("results/deconvolution/chunks/{chunk}/qc_within_sample.png"),
            # drivers_* deliberately NOT declared -> skipped by the script.
        params:
            samples        = SAMPLES,
            samples_subset = lambda w: _DECONV_CHUNKS[w.chunk],
            sample_bed_map = SAMPLE_FRAG_BED,
            genome_fasta   = config["genome"].get("fasta", ""),
            gene_annotation_bed = config["deconvolution"].get("gene_annotation_bed", ""),
            stage          = "within",
        conda:   "../../envs/python_analysis.yaml"
        threads: config["threads"]
        resources:
            mem_mb = _DECONV_WITHIN_MEM,
        script:  "../../scripts/module4_deconvolution.py"

    rule merge_within_sample:
        """Concatenate the per-chunk within-sample tables. Cheap."""
        input:
            scores  = expand("results/deconvolution/chunks/{chunk}/scores_within_sample.tsv",
                             chunk=_DECONV_CHUNK_IDS),
            qc      = expand("results/deconvolution/chunks/{chunk}/qc_within_sample.tsv",
                             chunk=_DECONV_CHUNK_IDS),
            windows = expand("results/deconvolution/chunks/{chunk}/background_signal_windows.bed.gz",
                             chunk=_DECONV_CHUNK_IDS),
        output:
            within_sample_scores      = "results/deconvolution/scores_within_sample.tsv",
            within_sample_qc          = "results/deconvolution/qc_within_sample.tsv",
            within_sample_windows_bed = "results/deconvolution/background_signal_windows.bed.gz",
        params:
            expected_samples = len(SAMPLES),
        conda:   "../../envs/python_analysis.yaml"
        threads: 1
        resources:
            mem_mb = 8000,
        script:  "../../scripts/module4_within_merge.py"
