"""
Module 7a — Batch × Group balance check.
Snakemake script: called via script: directive.

PURPOSE
Before running any statistical model, verify that the study groups
(Ctrl, NEW, MS-Rituximab-Stable, MS-Rituximab-Progressive) are not
completely confounded with sequencing/library preparation batches.

WHY THIS MATTERS
If all samples from a given group were processed in a single batch, then
any group-level difference in the count data could be due to batch effects
rather than biology. RUVg (Module 7b) can partially correct this, but it
cannot separate a signal that is 100% correlated with batch. Detecting
severe confounding early prevents over-interpreting results that are
really driven by technical variation.

OUTPUTS
  - A batch × group contingency table saved as TSV.
  - Printed warnings for each problematic configuration:
      * A single batch containing samples from only one group.
      * A group whose samples all come from a single batch.
  These warnings should be reported as limitations in the manuscript if
  confounding is present; they do NOT stop the pipeline.
"""

import pandas as pd
from pathlib import Path

# ── Step 1: Load metadata and restrict to QC-passing samples ─────────────────
# Samples with qc_include == FALSE (e.g. the QC_Mix pool) are excluded so they
# do not distort the batch/group contingency counts.
meta = pd.read_csv(snakemake.input.meta, sep="\t", comment="#")
meta = meta[meta["qc_include"].astype(str).str.upper() == "TRUE"]

batch_col = "batch"
group_col = "group"

# ── Step 2: Build the contingency table ──────────────────────────────────────
# pd.crosstab counts how many samples fall in each (batch, group) cell.
# Ideally, every group should appear in multiple batches and every batch
# should contain samples from multiple groups (balanced design).
table = pd.crosstab(meta[batch_col], meta[group_col])

# ── Step 3: Save the table ───────────────────────────────────────────────────
Path(snakemake.output.table).parent.mkdir(parents=True, exist_ok=True)
table.to_csv(snakemake.output.table, sep="\t")

print("[Module 7a] Batch × Group table:")
print(table.to_string())

# ── Step 4: Emit warnings for confounded configurations ──────────────────────
groups  = meta[group_col].unique()
batches = meta[batch_col].unique()

if len(batches) == 1 or all(meta[batch_col] == "NA"):
    # All samples in a single batch (or batch column not yet filled in).
    # RUVg will still run using constitutive anchors as negative controls,
    # but it will correct for sample-level technical noise rather than
    # discrete batch effects.
    print("[Module 7a] NOTE: All samples in a single batch or batch not assigned. "
          "RUVg will still run using constitutive anchors as negative controls.")
else:
    # Check for each batch: does it contain only one group?
    for b in batches:
        b_groups = meta.loc[meta[batch_col] == b, group_col].unique()
        if len(b_groups) == 1:
            print(
                f"[Module 7a] WARNING: batch '{b}' contains only group "
                f"'{b_groups[0]}' — severe batch/group confounding. "
                "Statistical correction will be incomplete; treat results cautiously."
            )
    # Check for each group: are all its samples from a single batch?
    for g in groups:
        g_batches = meta.loc[meta[group_col] == g, batch_col].unique()
        if len(g_batches) == 1 and len(batches) > 1:
            print(
                f"[Module 7a] WARNING: group '{g}' is entirely contained "
                f"in batch '{g_batches[0]}'. Any batch effect specific to "
                "this batch cannot be distinguished from a group effect."
            )
