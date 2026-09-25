"""
Module 4 — gather step for the chunked within-sample deconvolution.

Concatenates the per-chunk outputs of `deconvolve_within_chunk` back into the
canonical files the rest of the pipeline expects.

Why plain concatenation is correct here: every column in scores_within_sample
is computed from ONE sample in isolation.  `_ws_build_bg_model` reads only that
sample's own counts, alpha_s is fit from that sample's own background tiles,
and the BH-FDR is applied within each sample.  The cohort-median scaling factor
is never touched by the within-sample half.  So the union of the chunks is the
same set of rows the single-job path produces, only in chunk order.

background_signal_windows.bed.gz is the signal+background window set.  It is
built from the marker/off-target BEDs and does not depend on which samples were
scored, so every chunk writes an identical file and we keep the first.
"""

import gzip
import shutil

import pandas as pd

_scores_in  = list(snakemake.input.scores)
_qc_in      = list(snakemake.input.qc)
_windows_in = list(snakemake.input.windows)

print(f"[merge] {len(_scores_in)} chunks")

scores = pd.concat([pd.read_csv(p, sep="\t") for p in _scores_in], ignore_index=True)
qc     = pd.concat([pd.read_csv(p, sep="\t") for p in _qc_in], ignore_index=True)

# A sample landing in two chunks would silently double its rows downstream.
_dupes = scores.duplicated(subset=["sample", "cell_type"]).sum()
if _dupes:
    raise SystemExit(
        f"[merge] {_dupes} duplicated (sample, cell_type) rows — chunks overlap")

_n_samples = scores["sample"].nunique()
_expected  = int(snakemake.params.expected_samples)
if _n_samples != _expected:
    raise SystemExit(
        f"[merge] got {_n_samples} distinct samples, expected {_expected}; "
        "a chunk is missing or the sample list changed mid-run")

# Cohort-level permuted-signature guard, deferred here from
# module4_deconvolution.py. It has to run on the merged table: the rule is "every
# sample in the COHORT is flagged", and a chunk cannot know that. Same condition
# the single-job path applies, so both paths abort on exactly the same input.
#
# The distribution is printed either way because this metric runs close to its
# own threshold in normal operation (default 0.2; an accepted reference run had
# mean 0.199 and flagged 31% of samples), so "how far over" matters far more than
# the boolean when deciding whether a failure is real bias or mis-calibration.
#
# perm_flag is RECOMPUTED here rather than summed from the chunk column. The
# chunk column is written from the same config key, so on a clean run the two
# agree — but snakemake reuses up-to-date chunk outputs across a config change,
# so after retuning deconv_perm_thr only the merge reruns and the column still
# carries the OLD threshold's verdict. That is not hypothetical: 2026-09-21 the
# merge reported "550/550 flagged at threshold 0.6" while printing max=0.446,
# because the flags were made at 0.2 the day before. Deriving the boolean from
# median_perm_log2FE makes the message and the abort self-consistent, and means
# retuning the threshold does not require rerunning the deconvolution chunks.
_PERM_THR = float(snakemake.config.get("deconv_perm_thr", 0.2))

if "median_perm_log2FE" in qc.columns:
    _med = qc["median_perm_log2FE"]
    # Overwrite so the written QC table matches the threshold actually in force,
    # and so anything downstream filtering on perm_flag applies the same rule.
    # NaN > thr is False, so a sample with no permutation result counts as
    # unflagged — the same direction the chunk-written booleans took.
    qc["perm_flag"] = _med.abs() > _PERM_THR
    _n_flag = int(qc["perm_flag"].sum())
    print(f"[merge] permuted-signature null: {_n_flag}/{len(qc)} flagged "
          f"({_n_flag / len(qc):.1%}) at threshold {_PERM_THR:g}; "
          f"|median_perm_log2FE| "
          f"min={_med.abs().min():.3f} med={_med.abs().median():.3f} "
          f"p90={_med.abs().quantile(0.9):.3f} max={_med.abs().max():.3f}")
    if _n_flag == len(qc) and len(qc) > 0:
        raise SystemExit(
            f"[merge] ABORT: permuted-signature null fails for all {len(qc)} "
            f"samples (|median_perm_log2FE| median "
            f"{_med.abs().median():.3f} vs threshold {_PERM_THR:g}). Either a "
            "systematic GC/length bias (check RRE/DAC filters and genome FASTA) "
            "or a threshold set below this panel's intrinsic offset — the "
            "distribution above tells you which. Retune via the deconv_perm_thr "
            "config key; do NOT edit the default.")

scores.to_csv(snakemake.output.within_sample_scores, sep="\t", index=False)
qc.to_csv(snakemake.output.within_sample_qc, sep="\t", index=False)

# Identical across chunks by construction; verify rather than assume, since a
# mismatch would mean the chunks scored different window sets.
def _digest(path):
    with gzip.open(path, "rb") as fh:
        return hash(fh.read())

_ref = _digest(_windows_in[0])
for _p in _windows_in[1:]:
    if _digest(_p) != _ref:
        raise SystemExit(f"[merge] window set differs between chunks: {_p}")

shutil.copyfile(_windows_in[0], snakemake.output.within_sample_windows_bed)

print(f"[merge] {len(scores):,} rows / {_n_samples} samples "
      f"→ {snakemake.output.within_sample_scores}")
