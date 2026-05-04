"""
Module 6 — Assemble bin count matrix and apply mean-count filter.
Snakemake script: called via script: directive.

PURPOSE
Aggregate per-sample midpoint counts at the genome-tiling 5 kb bins
(produced by module6_make_bins.py) into a single matrix, then filter
out lowly-covered bins.

WHY FILTER BY MEAN COUNT?
Bins with very few counts across all samples carry essentially no signal
and only inflate the multiple-testing burden in DESeq2. Requiring a mean
count ≥ 5 across samples removes bins that are either in inaccessible
chromatin or outside the sequencing coverage. The threshold of 5 follows
standard practice (e.g. DESeq2's independent filtering step operates in
this range) and ensures each bin has enough counts for reliable
variance estimation.

INPUTS (via Snakemake)
  snakemake.input.counts      — ordered list of per-sample .counts files
  snakemake.input.bins        — the bin BED (output of module6_make_bins.py)
  snakemake.params.samples    — sample IDs in same order as counts
  snakemake.params.min_mean_count — minimum mean count threshold (default 5)

OUTPUT
  A TSV file (regions × samples) with only bins passing the mean-count filter.
  Used as the bin count matrix for DESeq2 and glmmTMB in Module 7.
"""

import pandas as pd
import numpy as np
from pathlib import Path

samples     = snakemake.params.samples
count_files = snakemake.input.counts
bins_bed    = snakemake.input.bins
min_mean    = float(snakemake.params.min_mean_count)

# ── Step 1: Define region IDs from the bin BED ───────────────────────────────
# Rows are identified as chr:start-end strings. This must match the order of
# rows in the per-sample .counts files, which were produced by running
# count_midpoints.py with the same bin BED file.
bins = pd.read_csv(bins_bed, sep="\t", header=None, usecols=[0, 1, 2])
bins.columns = ["chr", "start", "end"]
bins["region_id"] = (
    bins["chr"] + ":" +
    bins["start"].astype(str) + "-" +
    bins["end"].astype(str)
)

# ── Step 2: Stack per-sample count columns ───────────────────────────────────
# Identical to Module 5: the last column from each bedtools intersect -c file
# holds the integer count for that sample.
matrix = pd.DataFrame(index=bins["region_id"])
for sample, fpath in zip(samples, count_files):
    raw = pd.read_csv(fpath, sep="\t", header=None)
    matrix[sample] = raw.iloc[:, -1].values

# ── Step 3: Apply the mean-count filter ──────────────────────────────────────
# Compute the row-wise mean across all samples and keep only rows where this
# mean meets the threshold. Discarded rows will never appear in any downstream
# analysis (DESeq2, glmmTMB, ML).
row_means = matrix.mean(axis=1)
matrix    = matrix[row_means >= min_mean]

# ── Step 4: Write filtered matrix ────────────────────────────────────────────
Path(snakemake.output.matrix).parent.mkdir(parents=True, exist_ok=True)
matrix.to_csv(snakemake.output.matrix, sep="\t")

print(
    f"[Module 6] {matrix.shape[0]} bins (mean ≥ {min_mean}) "
    f"× {matrix.shape[1]} samples"
)
