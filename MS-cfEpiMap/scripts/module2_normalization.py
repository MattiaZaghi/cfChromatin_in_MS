"""
Module 2 — Compute constitutive anchor scaling factors.
Snakemake script: called via script: directive.

Inputs:  snakemake.input.counts  (list of per-sample .counts files)
         snakemake.input.anchors (constitutive anchors BED)
         snakemake.input.meta    (sample metadata TSV)
Outputs: snakemake.output.sf     (scaling factors TSV)
         snakemake.output.matrix (anchor count matrix TSV)

WHY THIS NORMALISATION?
DESeq2 computes its own size factors internally from the full count matrix, so we
do NOT pre-scale the counts going into differential analysis. However, for
everything outside DESeq2 — bigWig generation, B Cell Index computation,
deconvolution z-scores, TF AUC — we need a consistent per-sample scale factor
that removes sequencing-depth variation without distorting biology.

Constitutive anchor regions (Module 1 output) are H3K27ac peaks active in ≥14/16
reference cell types and are therefore expected to be equally active in every
sample. Any difference in counts at these regions between samples must be due to
library size differences, not biology. The scaling factor is simply the ratio of
the median anchor count (across all samples) to the per-sample anchor count:
    SF_i = median(anchor_total) / anchor_total_i
Multiplying sample i's raw counts by SF_i brings all samples to the same
effective library size.

The raw anchor count matrix produced here also serves as the negative-control
feature matrix for RUVg batch correction in Module 7b.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Step 1: Gather per-sample count files ─────────────────────────────────────
# snakemake.params.samples is the ordered list of sample IDs.
# snakemake.input.counts is the matching list of .counts files produced by
# count_midpoints.py with --regions constitutive_anchors.bed.
# We zip them to build a dict {sample_id: file_path} for ordered access.
samples     = snakemake.params.samples
count_files = {s: f for s, f in zip(samples, snakemake.input.counts)}

# ── Step 2: Read the constitutive anchor BED to define region IDs ─────────────
# Each anchor region gets a string ID of the form chr:start-end.
# This ID is used as the row index of the count matrix so downstream code
# can unambiguously identify which genomic interval each row corresponds to.
anchor_df = pd.read_csv(
    snakemake.input.anchors, sep="\t", header=None,
    usecols=[0, 1, 2],
    names=["chr", "start", "end"],
)
anchor_df["region_id"] = (
    anchor_df["chr"] + ":" +
    anchor_df["start"].astype(str) + "-" +
    anchor_df["end"].astype(str)
)
region_ids = anchor_df["region_id"].tolist()

# ── Step 3: Build the raw count matrix (regions × samples) ───────────────────
# Each .counts file was produced by bedtools intersect -c and has the same
# number of rows as the anchor BED, in the same order. The last column is the
# integer midpoint count. We stack one column per sample.
count_matrix = pd.DataFrame(index=region_ids)
for sample, fpath in count_files.items():
    raw = pd.read_csv(fpath, sep="\t", header=None)
    # The last column (iloc[:, -1]) is always the count from bedtools intersect -c
    count_matrix[sample] = raw.iloc[:, -1].values

# ── Step 4: Compute per-sample scaling factors ────────────────────────────────
# anchor_totals: sum of all anchor-region counts per sample (a proxy for total
# sequencing depth at constitutive loci).
# median_total: the target depth we scale every sample towards.
# scaling_factors: multipliers to be applied to raw counts in downstream tools.
#   A sample with twice the median depth gets SF = 0.5 (counts are halved).
#   A sample with half the median depth gets SF = 2.0 (counts are doubled).
anchor_totals   = count_matrix.sum(axis=0)
median_total    = anchor_totals.median()
scaling_factors = median_total / anchor_totals

# ── Step 5: Flag outlier samples ──────────────────────────────────────────────
# A scaling factor more than flag_sd standard deviations from the mean suggests
# the sample has very unusual library depth — either a failed library
# (very few anchor reads → very large SF) or contamination (many reads → tiny SF).
# These samples are NOT automatically excluded; they are flagged so the analyst
# can decide whether to investigate or remove them.
flag_sd   = float(snakemake.params.flag_sd)
sf_mean   = scaling_factors.mean()
sf_std    = scaling_factors.std()
flag_mask = (scaling_factors - sf_mean).abs() > flag_sd * sf_std

# ── Step 6: Load metadata and merge for group annotation in the output TSV ────
# Adding group and sample_type columns makes the SF table self-explanatory and
# directly useful for the QC barplot in generate_report.py.
meta = pd.read_csv(snakemake.input.meta, sep="\t", comment="#")

sf_df = pd.DataFrame({
    "sample_id":       anchor_totals.index,
    "anchor_reads":    anchor_totals.values,      # raw anchor count per sample
    "constitutive_sf": scaling_factors.values,    # multiplier for downstream normalisation
    "sf_flagged":      flag_mask.values,          # True = investigate this sample
})
sf_df = sf_df.merge(
    meta[["sample_id", "group", "sample_type", "protocol"]],
    on="sample_id", how="left",
)

# Print a warning for any flagged samples so they appear in the Snakemake log
if flag_mask.any():
    flagged = sf_df.loc[sf_df["sf_flagged"], "sample_id"].tolist()
    print(f"[Module 2] WARNING — {len(flagged)} samples flagged (SF > {flag_sd} SD): {flagged}")

# ── Step 7: Write outputs ─────────────────────────────────────────────────────
# sf output: TSV with one row per sample — used by Modules 3, 4, 8 for
#             constitutive-anchor-normalised density computations.
# matrix output: TSV with regions as rows and samples as columns — this raw
#                integer matrix is the negative-control input for RUVg in
#                Module 7b. It must NOT be pre-scaled.
Path(snakemake.output.sf).parent.mkdir(parents=True, exist_ok=True)
sf_df.to_csv(snakemake.output.sf, sep="\t", index=False)

Path(snakemake.output.matrix).parent.mkdir(parents=True, exist_ok=True)
count_matrix.to_csv(snakemake.output.matrix, sep="\t")

print(
    f"[Module 2] {len(samples)} samples | "
    f"median anchor reads: {int(median_total)} | "
    f"SF range: {scaling_factors.min():.3f} – {scaling_factors.max():.3f}"
)
