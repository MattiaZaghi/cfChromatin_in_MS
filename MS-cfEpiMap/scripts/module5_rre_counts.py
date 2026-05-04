"""
Module 5 — Assemble RRE count matrix.
Snakemake script: called via script: directive.

PURPOSE
Aggregate the per-sample midpoint count files (produced individually by
count_midpoints.py for each sample × region subset) into a single
regions × samples count matrix. This matrix is the primary input for
DESeq2 differential analysis in Module 7c.

WHY RAW COUNTS (NO PRE-SCALING)?
DESeq2 expects raw integer counts and derives its own size factors
internally using the median-ratio method. Pre-scaling the matrix would
distort the count distribution and invalidate DESeq2's statistical model.
The constitutive-anchor scaling factors from Module 2 are only used
outside DESeq2 (bigWigs, BCI, deconvolution).

INPUTS (via Snakemake)
  snakemake.input.counts   — ordered list of per-sample .counts files
  snakemake.input.regions  — BED file defining the RRE subset (universe,
                             CNS, immune, GWAS-proximal, or bin BED)
  snakemake.params.samples — list of sample IDs in the same order as counts
  snakemake.params.subset  — label for logging (e.g. "full", "cns", "immune")

OUTPUT
  A TSV file with region_id (chr:start-end) as row index and sample IDs as
  column headers, containing raw integer midpoint counts.
"""

import pandas as pd
from pathlib import Path

samples     = snakemake.params.samples
count_files = snakemake.input.counts   # list, same order as samples
regions_bed = snakemake.input.regions

# ── Step 1: Build the region ID index from the BED file ──────────────────────
# Each region gets a string ID of the form chr:start-end.
# The .counts files produced by count_midpoints.py have the same row order
# as the BED file they were run against, so we use this index directly.
regions = pd.read_csv(regions_bed, sep="\t", header=None, usecols=[0, 1, 2])
regions.columns = ["chr", "start", "end"]
regions["region_id"] = (
    regions["chr"] + ":" +
    regions["start"].astype(str) + "-" +
    regions["end"].astype(str)
)

# ── Step 2: Stack per-sample count columns into the matrix ───────────────────
# Each .counts file has the same number of rows as the BED (guaranteed by
# bedtools intersect -c, which outputs one line per -a record). The last
# column (iloc[:, -1]) is always the integer count appended by bedtools.
matrix = pd.DataFrame(index=regions["region_id"])
for sample, fpath in zip(samples, count_files):
    raw = pd.read_csv(fpath, sep="\t", header=None)
    matrix[sample] = raw.iloc[:, -1].values

# ── Step 3: Write the count matrix ───────────────────────────────────────────
# The output TSV is indexed by region_id (rows) with sample IDs as column names.
Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
matrix.to_csv(snakemake.output[0], sep="\t")

print(
    f"[Module 5] {snakemake.params.subset} RRE matrix: "
    f"{matrix.shape[0]} regions × {matrix.shape[1]} samples"
)
