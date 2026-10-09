configfile: "config.yaml"

import csv
import os
from pathlib import Path

# ─── Sample list ──────────────────────────────────────────────────────────────
DEFAULT_MARK = 'H3K27ac'

SAMPLES = []
SAMPLE_MARK = {}
if config.get('samplesheet'):
    samplesheet_path = Path(config['samplesheet'])
    if samplesheet_path.exists():
        with open(samplesheet_path, newline='') as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sid = row.get('sampleId') or row.get('sample') or row.get('sample_id')
                if sid:
                    SAMPLES.append(sid)
                    SAMPLE_MARK[sid] = (row.get('enrichment_mark') or DEFAULT_MARK).strip()
    else:
        print(f"Warning: samplesheet {samplesheet_path} not found; SAMPLES empty")
else:
    print("No samplesheet provided in config.yaml; define 'samplesheet'")

# Expose sample list to rules that need it (e.g. fragle_ct_estimation)
config['_samples_'] = SAMPLES
# Expose the per-sample enrichment mark so the enrichment score and QC rules
# can pick the matching chromHMM reference BEDs and cut-offs.
config['_sample_mark_'] = SAMPLE_MARK


def enrichment_refs(mark):
    """(ontarget, offtarget) chromHMM BEDs for one enrichment mark.

    Falls back to the single global chromhmm_ontarget/offtarget pair for any
    mark not listed under enrichment_targets, which keeps pre-existing configs
    (and the H3K27ac-only external cohort) working unchanged.
    """
    entry = config.get('enrichment_targets', {}).get(mark, {})
    return (
        entry.get('ontarget',
                  config.get('chromhmm_ontarget', 'ref_files/chromhmm_ontarget_H3K27ac.bed')),
        entry.get('offtarget',
                  config.get('chromhmm_offtarget', 'ref_files/chromhmm_offtarget_H3K27ac.bed')),
    )


def sample_mark(sample):
    return config.get('_sample_mark_', {}).get(sample, DEFAULT_MARK)

# ─── Include rule modules ──────────────────────────────────────────────────────
# Preprocessing
include: "rules/modules/fastqc.smk"
include: "rules/modules/trim.smk"

# Alignment
include: "rules/modules/align.smk"

# BAM processing chain
include: "rules/modules/sort_bam.smk"             # sort_bam, sort_readname_bam
include: "rules/modules/filter_properly_paired.smk"
include: "rules/modules/unique_sam.smk"
include: "rules/modules/quality_filter.smk"
include: "rules/modules/dedup.smk"
include: "rules/modules/dac_exclusion.smk"
include: "rules/modules/index_and_quality.smk"     # index_sam
include: "rules/modules/lib_complex_preseq.smk"
include: "rules/modules/createStatsSamtoolsfiltered.smk"
include: "rules/modules/snp_smash_fingerprint.smk"

# Fragments processing
include: "rules/modules/end_motif_gc.smk"
include: "rules/modules/calcFragsLength.smk"
include: "rules/modules/frag_size_plot.smk"
include: "rules/modules/bam_to_bed.smk"
include: "rules/modules/unique_frags.smk"
include: "rules/modules/frags_report.smk"

# Fragle CT estimation (optional)
include: "rules/modules/filter_bam_fragle.smk"
include: "rules/modules/fragle_ct_estimation.smk"
include: "rules/modules/ct_report.smk"

# Signal processing
include: "rules/modules/bedgraph_to_bigwig.smk"
include: "rules/modules/call_peaks.smk"
include: "rules/modules/peaks_annotations.smk"
include: "rules/modules/peaks_report.smk"
include: "rules/modules/enrichment.smk"
include: "rules/modules/enrichmentReport.smk"
include: "rules/modules/merge_enrichment_reports.smk"
include: "rules/modules/chromatin_count_normalization.smk"
include: "rules/modules/igv_reports.smk"

# Enrichment score QC (requires pre-built chromHMM reference BEDs)
include: "rules/modules/enrichment_score.smk"

# Sample QC flagging (data-driven knee-point thresholds)
include: "rules/modules/qc_flag.smk"

# External published cohort enrichment (Sadeh et al. GSM BED files)
include: "rules/modules/external_enrichment.smk"

# Meta plots
include: "rules/modules/meta_plot_housekeeping.smk"
include: "rules/modules/meta_plot_hk_promoters.smk"

# Reports
include: "rules/modules/quality_report_lite.smk"
include: "rules/modules/signal_report_lite.smk"
include: "rules/modules/merge_signal_reports.smk"
include: "rules/modules/multiqc.smk"

# Reference/download rules (used when ref files need to be fetched)
include: "rules/modules/download.smk"
include: "rules/modules/createGenomeIndex.smk"
include: "rules/modules/fetch_chrom_sizes.smk"

# ─── Rule ordering ────────────────────────────────────────────────────────────
ruleorder: trim_fastp > trim_galore

# ─── Helper: final processed BAM path (dac-filtered or dedup) ─────────────────
def get_final_bam(wildcards):
    base = config['outputFolder']
    if config.get('exclude_dac_regions', False):
        return f"{base}/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam"
    return f"{base}/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam"

def get_final_bai(wildcards):
    base = config['outputFolder']
    if config.get('exclude_dac_regions', False):
        return f"{base}/align/dac/{wildcards.sample}.dac_filtered.dedup.unique.sorted.bam.bai"
    return f"{base}/align/dedup/{wildcards.sample}.dedup.unique.sorted.bam.bai"

# ─── Targets ──────────────────────────────────────────────────────────────────
OUT = config['outputFolder']

DEDUP = expand(OUT + "/align/dedup/{sample}.dedup.unique.sorted.bam", sample=SAMPLES)

if config.get('exclude_dac_regions', False):
    FINAL_BAM = expand(OUT + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam", sample=SAMPLES)
    FINAL_BAI = expand(OUT + "/align/dac/{sample}.dac_filtered.dedup.unique.sorted.bam.bai", sample=SAMPLES)
else:
    FINAL_BAM = DEDUP
    FINAL_BAI = expand(OUT + "/align/dedup/{sample}.dedup.unique.sorted.bam.bai", sample=SAMPLES)

BED    = expand(OUT + "/frags/{sample}.bed",   sample=SAMPLES)
BIGWIG = expand(OUT + "/bigwig/{sample}.bw",   sample=SAMPLES)
PEAKS        = expand(OUT + "/peaks/{sample}.narrowPeak",        sample=SAMPLES)
FRAG_SIZES        = expand(OUT + "/frags/{sample}/{sample}.fragment_sizes.txt", sample=SAMPLES)
FRAG_SIZE_PERSAMPLE = expand(OUT + "/reports/qc/frag_size/{sample}_fragment_size.pdf", sample=SAMPLES)
FRAG_SIZE_GROUPS  = OUT + "/reports/qc/fragment_size_by_group.pdf"
FRAG_SIZE_PROTOCOL  = OUT + "/reports/qc/fragment_size_by_protocol.pdf"
MOTIFS       = expand(OUT + "/motifs/{sample}/{sample}_4NMER_bp_motif.bed", sample=SAMPLES)
LIB_COMPLEX  = expand(OUT + "/align/{sample}/{sample}.lc_extrap.txt", sample=SAMPLES)
UNIQUE_FRAGS = expand(OUT + "/frags/{sample}/{sample}_unique_frags.csv", sample=SAMPLES)

META_PLOTS    = expand(OUT + "/reports/meta_plots/{sample}_housekeeping_meta.pdf", sample=SAMPLES)
HK_PROM_META  = expand(OUT + "/reports/qc/hk_meta/{sample}_hk_promoters_meta.pdf", sample=SAMPLES)

# Enrichment score targets — require a usable reference pair for EVERY mark the
# samplesheet actually uses, so a missing H3K4me3 BED disables the block loudly
# instead of silently scoring H3K4me3 libraries against H3K27ac regions.
def _usable(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

_MARKS_IN_USE = sorted(set(SAMPLE_MARK.values())) or [DEFAULT_MARK]
_MISSING_REFS = {
    m: [p for p in enrichment_refs(m) if not _usable(p)]
    for m in _MARKS_IN_USE
}
_MISSING_REFS = {m: p for m, p in _MISSING_REFS.items() if p}
for _m, _paths in _MISSING_REFS.items():
    print(f"Warning: enrichment reference BEDs missing/empty for mark {_m}: "
          f"{', '.join(_paths)} — enrichment score targets disabled")
_ENRICHMENT_AVAILABLE = not _MISSING_REFS
ENRICHMENT_PLOTS        = expand(OUT + "/reports/enrichment/{sample}_enrichment_score.pdf", sample=SAMPLES)
ENRICHMENT_COMPARISON   = OUT + "/reports/enrichment/enrichment_comparison_plots.pdf"
ENRICHMENT_CROSS_COMP   = OUT + "/reports/enrichment/enrichment_cross_comparison.pdf"

MULTIQC_FASTQC = OUT + "/reports/multiqc/multiqc_fastqc.html"

TARGETS = []
TARGETS.extend(FINAL_BAM)
TARGETS.extend(FINAL_BAI)
TARGETS.extend(BED)
TARGETS.extend(BIGWIG)
TARGETS.extend(PEAKS)
TARGETS.extend(FRAG_SIZES)
TARGETS.extend(FRAG_SIZE_PERSAMPLE)
TARGETS.append(FRAG_SIZE_GROUPS)
TARGETS.append(FRAG_SIZE_PROTOCOL)
TARGETS.extend(UNIQUE_FRAGS)
TARGETS.extend(LIB_COMPLEX)
TARGETS.extend(META_PLOTS)
TARGETS.extend(HK_PROM_META)
TARGETS.append(MULTIQC_FASTQC)

if _ENRICHMENT_AVAILABLE and SAMPLES:
    TARGETS.extend(ENRICHMENT_PLOTS)
    TARGETS.append(ENRICHMENT_COMPARISON)
    TARGETS.append(OUT + "/qc/qc_summary.tsv")
    TARGETS.append(OUT + "/qc/qc_summary_protocol.tsv")
    TARGETS.append(OUT + "/reports/qc/qc_distributions.pdf")
    TARGETS.append(OUT + "/reports/qc/qc_distributions_protocol.pdf")
    # Cross-comparison with external GSM cohort (only if external BED dir exists)
    _EXT_BED_DIR = config.get('external_bed_dir', '/date/gcb/gcb_MZ/Analysis/BED/H3K27ac')
    if os.path.isdir(_EXT_BED_DIR):
        TARGETS.append(ENRICHMENT_CROSS_COMP)

if config.get('read_method', 'PE') == 'PE':
    TARGETS.extend(MOTIFS)

if config.get('fragle_ct_estimation', False) and SAMPLES:
    TARGETS.append(OUT + "/reports/fragle/Fragle.txt")
    TARGETS.append(OUT + "/reports/multiqc/ct_fragle_mqc.csv")

rule all:
    input: TARGETS
