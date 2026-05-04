"""
Module 3 — B Cell Index and Reconstitution Analysis.
Snakemake script: called via script: directive.

WHY A B CELL INDEX?
Rituximab is an anti-CD20 monoclonal antibody that depletes circulating B cells.
If B cell reconstitution after Rituximab varies between patients, the cfDNA pool
will contain different amounts of B-cell-derived chromatin fragments. Without
accounting for this, any group-level difference we see in H3K27ac signal at
B-cell-active regions could simply reflect reconstitution status rather than
disease biology.

The B Cell Index (BCI) quantifies B-cell-derived cfDNA per sample:
    BCI = (midpoints/kb at bcell_rre.bed) / (midpoints/kb at constitutive_anchors.bed)

The ratio form is important: it is internally self-normalised because numerator
and denominator are computed from the same fragment BED file, so there is no need
to apply the Module 2 scaling factor here.

After computing BCI, we:
  1. Scale it by dividing by the Ctrl group mean → bci_scaled (used as a covariate
     in the DESeq2 design formula in Module 7c).
  2. Run Mann-Whitney U tests between all group pairs to see if BCI differs.
  3. Correlate BCI with measured CD19+ B-cell counts (assay validation).
  4. Correlate BCI with months since last Rituximab infusion in treated patients
     to characterise the reconstitution time course.
  5. Flag treated samples with BCI > 1.5× Ctrl mean for review.
"""

import subprocess
import tempfile
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")   # non-interactive backend required on compute nodes
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import seaborn as sns


# ─── Utility functions ────────────────────────────────────────────────────────

def kb_coverage(bed_path: str) -> float:
    """
    Return the total kilobases spanned by all intervals in a BED file.
    Used to compute fragment-per-kilobase density so that regions of
    different sizes are compared on an equal footing.
    """
    total = 0
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                total += int(parts[2]) - int(parts[1])
    return total / 1000.0


def count_midpoints_in_bed(frags_bed: str, regions_bed: str) -> int:
    """
    Count the total number of fragment midpoints that fall within any
    interval in regions_bed.

    Steps:
      1. Convert each fragment (chr, start, end) → 1-bp midpoint (chr, mid, mid+1).
      2. Sort the midpoints for efficient bedtools sweep.
      3. Run bedtools intersect -c to count midpoints per region.
      4. Sum all per-region counts and return the total.

    The temporary file is created in /tmp and cleaned up on exit.
    """
    # Step 1+2: compute midpoints and sort via a shell pipeline
    awk = (
        "awk 'BEGIN{OFS=\"\\t\"} "
        "{mid=int(($2+$3)/2); print $1,mid,mid+1}' "
        f"{frags_bed}"
    )
    with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        p1 = subprocess.Popen(awk, shell=True, stdout=subprocess.PIPE)
        p2 = subprocess.run(
            "bedtools sort -i -", shell=True,
            stdin=p1.stdout, capture_output=True, text=False,
        )
        p1.stdout.close()
        with open(tmp_path, "wb") as f:
            f.write(p2.stdout)

        # Step 3: count midpoints per region, then sum across all regions (step 4)
        res = subprocess.run(
            f"bedtools intersect -a {regions_bed} -b {tmp_path} -c",
            shell=True, capture_output=True, text=True,
        )
        total = sum(int(l.split()[-1]) for l in res.stdout.splitlines() if l.strip())
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return total


# ─── Load Snakemake inputs and parameters ────────────────────────────────────
samples   = snakemake.params.samples      # list of sample IDs to process
frags_dir = snakemake.params.frags_dir   # directory containing {sample_id}.bed files
bcell_bed = snakemake.input.bcell_rre    # BED of B-cell-specific regulatory elements
anchors   = snakemake.input.anchors      # constitutive anchor BED (denominator)
flag_fold = float(snakemake.params.flag_fold)  # fold-threshold for flagging reconstitution

meta = pd.read_csv(snakemake.input.meta, sep="\t", comment="#")
meta = meta[meta["sample_id"].isin(samples)].copy()

# Pre-compute kb coverages (constant across samples — calculated once)
bcell_kb  = kb_coverage(bcell_bed)
anchor_kb = kb_coverage(anchors)

# ─── Step 1: Compute raw BCI for every sample ────────────────────────────────
# For each sample, count midpoints in the B-cell RRE set and in the constitutive
# anchor set, then form the density ratio.
records = []
for sid in samples:
    frags = f"{frags_dir}/{sid}.bed"
    b_cnt = count_midpoints_in_bed(frags, bcell_bed)    # B-cell-region midpoints
    a_cnt = count_midpoints_in_bed(frags, anchors)       # anchor-region midpoints
    # BCI = (B-cell density) / (anchor density)  [dimensionless ratio]
    bci   = (b_cnt / bcell_kb) / (a_cnt / anchor_kb) if a_cnt > 0 else np.nan
    records.append({"sample_id": sid, "bcell_counts": b_cnt,
                    "anchor_counts": a_cnt, "bci": bci})
    print(f"  {sid}: BCI={bci:.4f}")

# ─── Step 2: Scale BCI relative to the Ctrl group mean ───────────────────────
# bci_scaled = BCI / mean(BCI in Ctrl samples)
# This makes the Ctrl group average = 1.0 and treated patients are expressed
# as multiples of healthy baseline. This scaled value enters the DESeq2 design
# as a continuous covariate (Module 7c) to partial out B-cell reconstitution
# variance before testing for group differences at other regulatory elements.
bci_df = pd.DataFrame(records)
bci_df = bci_df.merge(meta[["sample_id", "group", "sample_type",
                             "cd19_count", "months_since_rtx"]], on="sample_id")

ctrl_mean = bci_df.loc[bci_df["group"] == "Ctrl", "bci"].mean()
bci_df["bci_scaled"] = bci_df["bci"] / ctrl_mean

# ─── Step 3: Flag treated samples with unexpectedly high B cell reconstitution ─
# A treated sample with BCI > flag_fold × Ctrl mean has reconstituted B cells
# beyond the expected range. This is not necessarily exclusion-worthy — rapid
# reconstitution may itself be biologically meaningful — but it warrants
# manual cross-check with CD19+ count and months_since_rtx.
treated = bci_df["group"].isin(["MS-Rituximab-Stable", "MS-Rituximab-Progressive"])
bci_df["bci_flag"] = treated & (bci_df["bci"] > flag_fold * ctrl_mean)

# ─── Step 4: Pairwise Mann-Whitney U tests between groups ────────────────────
# Non-parametric test: does BCI differ significantly between any pair of groups?
# The MS-Rituximab-Stable vs MS-Rituximab-Progressive comparison is the critical
# one: higher BCI in Progressive patients would suggest that faster B cell
# reconstitution is associated with non-response to Rituximab.
groups    = ["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
stat_rows = []
for i, g1 in enumerate(groups):
    for g2 in groups[i + 1:]:
        a = bci_df.loc[bci_df["group"] == g1, "bci"].dropna()
        b = bci_df.loc[bci_df["group"] == g2, "bci"].dropna()
        if len(a) > 1 and len(b) > 1:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            stat_rows.append({"group1": g1, "group2": g2, "U": u, "pvalue": p})
# Store test results as a DataFrame attribute for potential downstream access
bci_df.attrs["mw_tests"] = pd.DataFrame(stat_rows)

# ─── Step 5: Spearman correlation — BCI vs measured CD19+ count ───────────────
# This validates the BCI assay: if cfDNA B-cell signal truly reflects circulating
# B cells, BCI should correlate strongly with the direct flow-cytometry CD19+ count
# (expected r > 0.5 based on Sadeh et al. original cfChIP-seq validation).
cd19_df = bci_df[["bci", "cd19_count"]].replace("NA", np.nan).dropna()
if len(cd19_df) > 5:
    r_cd19, p_cd19 = stats.spearmanr(cd19_df["bci"], cd19_df["cd19_count"].astype(float))
    print(f"[Module 3] BCI vs CD19+: r={r_cd19:.3f}, p={p_cd19:.4f}")

# ─── Step 6: Spearman correlation — BCI vs months since Rituximab ─────────────
# In treated patients, B cells reconstitute over months after each Rituximab
# infusion. Correlating BCI with time since last infusion quantifies this
# reconstitution trajectory from cfDNA rather than blood counts.
rtx_df = bci_df[bci_df["group"].isin(
    ["MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
)][["bci", "months_since_rtx"]].replace("NA", np.nan).dropna()
if len(rtx_df) > 5:
    r_rtx, p_rtx = stats.spearmanr(rtx_df["bci"], rtx_df["months_since_rtx"].astype(float))
    print(f"[Module 3] BCI vs months_since_rtx: r={r_rtx:.3f}, p={p_rtx:.4f}")

# ─── Step 7: Save outputs ────────────────────────────────────────────────────
Path(snakemake.output.scores).parent.mkdir(parents=True, exist_ok=True)
bci_df.to_csv(snakemake.output.scores, sep="\t", index=False)

# Write a separate file listing only flagged samples for the HTML report
flagged = bci_df[bci_df["bci_flag"]]
flagged.to_csv(snakemake.output.flagged, sep="\t", index=False)
if not flagged.empty:
    print(f"[Module 3] {len(flagged)} treated samples flagged BCI > {flag_fold}× Ctrl mean")

# ─── Step 8: Figures ─────────────────────────────────────────────────────────
# Consistent colour palette used across all modules
palette = {
    "Ctrl":                      "#4878CF",
    "NEW":                       "#6ACC65",
    "MS-Rituximab-Stable":       "#D65F5F",
    "MS-Rituximab-Progressive":  "#B47CC7",
}

with pdf_backend.PdfPages(snakemake.output.plots) as pp:

    # ── Figure 1: Violin plot — BCI distribution by group ────────────────────
    # The violin shows the full distribution shape; the stripplot overlays
    # individual data points so that small group sizes do not hide outliers.
    fig, ax = plt.subplots(figsize=(7, 5))
    order = [g for g in groups if g in bci_df["group"].unique()]
    sns.violinplot(data=bci_df, x="group", y="bci", order=order,
                   palette=palette, ax=ax, inner="box", cut=0)
    sns.stripplot(data=bci_df, x="group", y="bci", order=order,
                  palette=palette, ax=ax, size=4, alpha=0.7, jitter=True)
    ax.set_title("B Cell Index by group")
    ax.set_xlabel("")
    ax.set_ylabel("BCI")
    ax.tick_params(axis="x", rotation=30)
    plt.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

    # ── Figure 2: BCI vs CD19+ scatter — assay validation ───────────────────
    # Each point is a sample; colour = group. A tight correlation validates
    # that the cfDNA B-cell signal tracks peripheral B cell counts.
    if len(cd19_df) > 5:
        fig, ax = plt.subplots(figsize=(5, 5))
        merged = bci_df[["sample_id", "group", "bci", "cd19_count"]].copy()
        merged["cd19_count"] = pd.to_numeric(merged["cd19_count"], errors="coerce")
        merged = merged.dropna(subset=["bci", "cd19_count"])
        for grp, sub in merged.groupby("group"):
            ax.scatter(sub["cd19_count"], sub["bci"],
                       label=grp, color=palette.get(grp, "grey"), alpha=0.8)
        ax.set_xlabel("CD19+ count (cells/µl)")
        ax.set_ylabel("BCI")
        ax.set_title(f"BCI vs CD19+  (r={r_cd19:.2f}, p={p_cd19:.3f})")
        ax.legend(fontsize=7)
        plt.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

    # ── Figure 3: BCI vs months since Rituximab — reconstitution time course ─
    # Only Rituximab-treated patients are plotted. An upward trend over time
    # would confirm that BCI captures B-cell reconstitution dynamics in cfDNA.
    if len(rtx_df) > 5:
        fig, ax = plt.subplots(figsize=(5, 5))
        tr = bci_df[bci_df["group"].isin(
            ["MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
        )].copy()
        tr["months_since_rtx"] = pd.to_numeric(tr["months_since_rtx"], errors="coerce")
        tr = tr.dropna(subset=["bci", "months_since_rtx"])
        for grp, sub in tr.groupby("group"):
            ax.scatter(sub["months_since_rtx"], sub["bci"],
                       label=grp, color=palette.get(grp, "grey"), alpha=0.8)
        ax.set_xlabel("Months since last Rituximab infusion")
        ax.set_ylabel("BCI")
        ax.set_title(f"B cell reconstitution  (r={r_rtx:.2f}, p={p_rtx:.3f})")
        ax.legend(fontsize=7)
        plt.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

print("[Module 3] Done.")
