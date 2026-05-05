"""
Module 4 — Cell-Type Deconvolution (optimized).

KEY OPTIMIZATION vs previous version:
  Old: 16 cell types × 56 samples = 896 serial bedtools calls, each sorting
       the full fragment BED from scratch.
  New: Build one combined regions BED (all cell types labeled). Sort each
       sample's midpoints ONCE, then one bedtools intersect covers all cell
       types. Parallelized across samples → ~56 parallel jobs total.
       Typical runtime: 5-15 min instead of 6+ hours.
"""

import os
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import seaborn as sns


# ─── Cell-type definitions ────────────────────────────────────────────────────
CELL_TYPES = {
    "Oligodendrocyte": "cns",
    "OPC":             "cns",
    "Neuron":          "cns",
    "Astrocyte":       "cns",
    "Microglia":       "cns",
    "CD4_Th1":         "immune",
    "CD4_Th17":        "immune",
    "CD4_Treg":        "immune",
    "CD8_T":           "immune",
    "B_naive":         "bcell",
    "B_memory":        "bcell",
    "Monocyte_CD14":   "immune",
    "Monocyte_CD16":   "immune",
    "NK":              "immune",
    "Neutrophil":      "immune",
    "Megakaryocyte":   "other",
}

CDI_TYPES = ["Oligodendrocyte", "OPC", "Neuron"]
NII_TYPES = ["CD4_Th1", "CD4_Th17", "Monocyte_CD14", "NK"]
IAI_TYPES = ["Monocyte_CD14", "Monocyte_CD16", "NK"]


# ─── Parallel worker (module-level for pickling) ──────────────────────────────

def _count_one_sample(args):
    """
    Sort midpoints ONCE for one sample, then intersect with the combined
    regions BED (all cell types in one file, labeled in col 4).
    Returns (sample_id, {cell_type: total_count}).
    """
    sid, frags_bed, combined_bed = args
    mid_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    try:
        awk = (
            f"awk 'BEGIN{{OFS=\"\\t\"}}"
            f"{{mid=int(($2+$3)/2); print $1,mid,mid+1}}' {frags_bed}"
        )
        p1 = subprocess.Popen(awk, shell=True, stdout=subprocess.PIPE)
        p2 = subprocess.run(
            "bedtools sort -i -", shell=True,
            stdin=p1.stdout, capture_output=True,
        )
        p1.stdout.close()
        with open(mid_tmp, "wb") as f:
            f.write(p2.stdout)

        res = subprocess.run(
            f"bedtools intersect -a {combined_bed} -b {mid_tmp} -c",
            shell=True, capture_output=True, text=True,
        )
    finally:
        if os.path.exists(mid_tmp):
            os.remove(mid_tmp)

    ct_totals: dict[str, int] = defaultdict(int)
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 5:
            ct_totals[parts[3]] += int(parts[-1])

    return sid, dict(ct_totals)


# ─── Step 1: Load inputs ──────────────────────────────────────────────────────
samples   = snakemake.params.samples
frags_dir = snakemake.params.frags_dir
n_cpus    = snakemake.threads

meta  = pd.read_csv(snakemake.input.meta, sep="\t", comment="#")
meta  = meta[meta["sample_id"].isin(samples)]
sf_df = pd.read_csv(snakemake.input.sf, sep="\t")
sf    = sf_df.set_index("sample_id")["constitutive_sf"].to_dict()

# anchor_counts: raw fragment midpoint counts at constitutive anchor regions per
# sample, produced by Module 2 alongside the scaling factors.  Used as the
# within-sample null denominator for the Poisson enrichment test.
anchor_counts = sf_df.set_index("sample_id")["anchor_reads"].to_dict()

# anchor_kb: total kilobases spanned by all constitutive anchor regions.
# Computed once here so it does not need to be re-read inside the per-cell-type loop.
# We open the same BED file that was used by Module 2 so the denominator is
# guaranteed to be consistent with anchor_counts above.
anchor_kb = sum(
    int(line.split()[2]) - int(line.split()[1])
    for line in open(snakemake.input.anchors)
    if not line.startswith("#") and len(line.split()) >= 3
) / 1000.0

# ─── Step 1b: Estimate negative-binomial dispersion from anchor counts ────────
#
# WHY NEGATIVE BINOMIAL INSTEAD OF POISSON?
# The Poisson distribution assumes variance = mean.  Real sequencing count data
# is almost always overdispersed: the observed variance is larger than the mean
# because of PCR amplification bias, GC-content effects, and non-uniform
# chromatin accessibility.  The negative-binomial (NB) distribution adds a
# single extra parameter α (dispersion) that allows variance > mean:
#
#   Var(X) = μ + α · μ²
#
# When α → 0 the NB collapses back to Poisson, so the Poisson is a special case.
#
# HOW WE ESTIMATE α FROM THE ANCHOR REGIONS
# The constitutive anchor regions are biologically invariant across all samples
# (they are active in ≥14/16 reference cell types by definition).  Any
# sample-to-sample variability in their counts is therefore purely technical:
# PCR stochasticity, GC-content sequencing bias, IP efficiency fluctuations.
# That technical variability is exactly what overdispersion models.
# By computing α from anchor counts we get an experiment-specific noise floor
# that is free of biological signal — precisely what we want for a null model.
#
# STEP-BY-STEP PROCEDURE
#
# 1. Load the raw anchor count matrix (regions × samples) from Module 2.
#    This matrix contains integer fragment midpoint counts at every constitutive
#    anchor region for every sample.  It was NOT pre-scaled — we need raw counts
#    here so that the depth normalisation below is mathematically correct.

anchor_mat = pd.read_csv(snakemake.input.anchor_matrix, sep="\t", index_col=0)

# Keep only the samples that are part of this run (the matrix may contain more).
anchor_mat = anchor_mat[[s for s in samples if s in anchor_mat.columns]]

# 2. Depth-normalise each sample's counts.
#
#    "Depth" here means the total number of fragment midpoints that landed in
#    anchor regions for that sample (i.e. the column sum).  We divide each
#    column by its sum and then multiply by the median column sum across all
#    samples.  After this operation every sample has the same total anchor count
#    (the median), so comparing variances across samples is fair — differences
#    no longer reflect sequencing depth but only biological/technical noise.
#
#    Concretely for sample s and region r:
#      normalised_count(r, s) = raw_count(r, s)
#                               ────────────────────────────── × median_depth
#                               total_anchor_count(s)
#
#    This is the same principle as CPM (counts per million) normalisation but
#    scaled to the median depth rather than 1 million so the counts stay in a
#    range similar to the raw values, which makes the dispersion estimate
#    directly comparable to the raw counts used in the enrichment test.

col_sums   = anchor_mat.sum(axis=0)           # total anchor reads per sample
median_depth = float(col_sums.median())        # target depth: median across cohort

# Avoid division by zero for any sample with 0 anchor reads (degenerate case).
col_sums_safe = col_sums.replace(0, np.nan)

# Broadcast division: divide each element by its column sum, then scale up.
# The result has the same shape as anchor_mat (regions × samples).
anchor_norm = anchor_mat.div(col_sums_safe, axis=1) * median_depth

# 3. For each anchor region compute the mean and variance across all samples.
#
#    We use axis=1 (across columns = across samples) so we get one number per
#    row (region).  ddof=1 gives the unbiased sample variance (divides by N-1
#    rather than N), which is correct when estimating population variance from
#    a finite sample.

region_mean = anchor_norm.mean(axis=1)         # shape: (n_regions,)
region_var  = anchor_norm.var(axis=1, ddof=1)  # shape: (n_regions,)

# 4. Derive the dispersion α for each region.
#
#    From the NB variance formula  Var = μ + α·μ²  we rearrange:
#      α = (Var − μ) / μ²
#
#    This can be negative when the observed variance is smaller than the mean
#    (i.e. the data are underdispersed relative to Poisson — rare but possible
#    in very clean experiments or with few samples).  We clip at a small
#    positive floor (1e-6) rather than zero so that the NB distribution is
#    always well-defined.  A floor of 1e-6 is functionally equivalent to
#    Poisson: variance ≈ mean + 1e-6 · mean² ≈ mean for realistic count ranges.
#
#    We also mask regions where mean ≈ 0 (mean < 1) because α is undefined when
#    the denominator μ² ≈ 0 — these are regions with near-zero signal that
#    add noise to the dispersion estimate.

valid_regions = region_mean >= 1.0             # boolean mask: skip near-zero regions
alpha_per_region = (
    (region_var[valid_regions] - region_mean[valid_regions])
    / (region_mean[valid_regions] ** 2)
)

# 5. Aggregate to a single global dispersion estimate.
#
#    We take the MEDIAN rather than the mean for robustness: a small number of
#    highly variable anchor regions (e.g. near repeat elements or ENCODE
#    blacklist-adjacent sites) could pull the mean far upward, leading to an
#    overly conservative test for every cell type.  The median is insensitive
#    to such outliers and represents the "typical" technical noise level.
#
#    The floor of 1e-6 is applied after aggregation as a safety net.

nb_alpha = float(np.median(alpha_per_region.clip(lower=0)))
nb_alpha = max(nb_alpha, 1e-6)   # ensure strict positivity

print(
    f"[Module 4] NB dispersion estimated from {valid_regions.sum()} anchor regions "
    f"(median α = {nb_alpha:.6f}).  "
    f"{'α ≈ 0: NB ≈ Poisson' if nb_alpha < 1e-4 else 'Overdispersion detected.'}"
)

# ─── Step 2: Load RRE universe and define signature regions ───────────────────
rre = pd.read_csv(snakemake.input.rre_universe, sep="\t", header=None)
rre.columns = list(range(rre.shape[1]))
rre["cell_types"] = rre[3].fillna("").astype(str) if rre.shape[1] >= 4 else ""

ct_regions: dict[str, pd.DataFrame] = {}
for ct in CELL_TYPES:
    in_ct    = rre["cell_types"].str.contains(ct, regex=False)
    others   = [c for c in CELL_TYPES if c != ct]
    in_other = rre["cell_types"].apply(lambda x: any(c in x for c in others))
    ct_regions[ct] = rre[in_ct & ~in_other]

# sig_kb: total kilobases of exclusive signature regions for each cell type.
# This is the target space size used in the Poisson null — how much genomic
# territory we are counting over for each cell type.
# rre columns 1 and 2 are start and end (integer); subtracting and summing
# gives total base pairs, divided by 1000 gives kilobases.
sig_kb: dict[str, float] = {
    ct: float((sig_df[2] - sig_df[1]).sum()) / 1000.0
    for ct, sig_df in ct_regions.items()
}

# ─── Step 3: Write combined BED (all cell types, label in col 4) ──────────────
combined_rows = []
for ct, sig_df in ct_regions.items():
    if sig_df.empty:
        print(f"  [Module 4] {ct}: 0 signature regions — skipping")
        continue
    print(f"  [Module 4] {ct}: {len(sig_df)} signature regions")
    for _, row in sig_df.iterrows():
        combined_rows.append(f"{row[0]}\t{int(row[1])}\t{int(row[2])}\t{ct}")

raw_combined = tempfile.NamedTemporaryFile(
    suffix=".bed", delete=False, mode="w"
)
raw_combined.write("\n".join(combined_rows) + "\n")
raw_combined.close()

sorted_combined = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
subprocess.run(
    f"bedtools sort -i {raw_combined.name} > {sorted_combined}",
    shell=True, check=True,
)
os.remove(raw_combined.name)

# ─── Step 4: Count midpoints per sample in parallel ───────────────────────────
print(f"\n[Module 4] Counting for {len(samples)} samples "
      f"using {n_cpus} parallel workers ...")

args_list = [
    (sid, f"{frags_dir}/{sid}.bed", sorted_combined)
    for sid in samples
]

sig_counts: dict[str, dict[str, int]] = defaultdict(dict)
with ProcessPoolExecutor(max_workers=n_cpus) as pool:
    for sid, ct_cnts in pool.map(_count_one_sample, args_list):
        for ct in CELL_TYPES:
            sig_counts[ct][sid] = ct_cnts.get(ct, 0)
        print(f"  done: {sid}")

os.remove(sorted_combined)

# ─── Step 5: Anchor-normalise, Poisson enrichment test, and z-scores ─────────
ctrl_samples = meta.loc[meta["group"] == "Ctrl", "sample_id"].tolist()


def _bh_correct(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction on a 1-D array of p-values.
    Returns adjusted p-values (padj) clipped to [0, 1].

    How it works step by step:
      1. Sort p-values from smallest to largest and remember the original
         positions so we can unsort at the end.
      2. Multiply each sorted p-value by n / rank (where rank goes 1..n).
         This is the BH formula: padj_i = p_i * n / rank_i.
      3. Enforce monotonicity from right to left with a cumulative minimum:
         a later (larger) p-value must never produce a smaller padj than an
         earlier one, otherwise the rejection region would not be contiguous.
      4. Put the adjusted values back in original order and clip to 1.
    """
    n = len(pvals)
    if n == 0:
        return np.array([])
    # Step 1: argsort gives positions that would sort the array
    order   = np.argsort(pvals)
    # Step 2: multiply by n / rank (ranks are 1-based: 1, 2, ..., n)
    bh      = pvals[order] * n / np.arange(1, n + 1)
    # Step 3: cumulative minimum from right — np.minimum.accumulate on reversed
    # array, then reverse back
    padj_sorted = np.minimum.accumulate(bh[::-1])[::-1]
    # Step 4: unsort back to original order
    padj = np.empty(n)
    padj[order] = padj_sorted
    return np.minimum(padj, 1.0)


rows = []
for ct, sample_cnts in sig_counts.items():

    # ── Group-level reference (same as before) ────────────────────────────────
    # Collect anchor-normalised counts from all control samples for this cell
    # type.  ctrl_mean is the reference level: the average normalised signal
    # healthy controls show at this cell type's signature.
    ctrl_norms = [
        sample_cnts[s] * sf.get(s, 1.0)
        for s in ctrl_samples if s in sample_cnts
    ]
    ctrl_mean = np.mean(ctrl_norms) if ctrl_norms else 1.0

    # ── Per-sample Poisson null parameters ───────────────────────────────────
    # s_kb: total kilobases of exclusive signature regions for this cell type.
    # Used to scale the anchor density to get the expected count.
    s_kb = sig_kb.get(ct, 0.0)

    for sid, cnt in sample_cnts.items():

        # ── Anchor-normalised count and group-relative z-score ────────────────
        # norm_cnt: raw count scaled by this sample's constitutive anchor
        # scaling factor, placing it on a common library-size-adjusted scale.
        norm_cnt = cnt * sf.get(sid, 1.0)

        # z: how many Poisson standard deviations this sample sits above or
        # below the control mean.  sqrt(ctrl_mean) is the Poisson SD.
        z = (norm_cnt - ctrl_mean) / np.sqrt(max(ctrl_mean, 1))

        # ── Within-sample Poisson enrichment test ────────────────────────────
        # anchor_cnt: raw fragment midpoints at constitutive anchor regions for
        # this specific sample (from Module 2 output, not normalised).
        anchor_cnt = anchor_counts.get(sid, 0)

        if anchor_cnt > 0 and anchor_kb > 0 and s_kb > 0:

            # anchor_density: reads per kilobase at constitutive loci in this
            # specific sample.  Both numerator and denominator are raw counts
            # from the same library, so this ratio is automatically normalised
            # for sequencing depth — a deeper library has proportionally more
            # reads in both numerator (anchor_cnt) and denominator (anchor_kb
            # is fixed), so the density tracks IP efficiency rather than depth.
            anchor_density = anchor_cnt / anchor_kb

            # expected: the null-hypothesis count at this cell type's signature
            # regions.  If cfDNA fragments were distributed at the same density
            # as the constitutive anchors (no cell-type-specific enrichment),
            # you would observe this many midpoints.
            expected = anchor_density * s_kb

            # fold_enrichment: how many times above the null expectation the
            # observed count is.  A value of 1 means no enrichment; 2 means
            # twice the constitutive baseline; <1 means depletion.
            fold_enrichment = cnt / expected if expected > 0 else np.nan

            # ── Negative-binomial one-sided enrichment test ───────────────────
            #
            # We now use the NB distribution rather than Poisson.  The key
            # difference is that NB has an extra parameter α (nb_alpha,
            # estimated from anchor regions in Step 1b) that allows the
            # variance to exceed the mean:
            #
            #   Var(X_NB) = μ + α · μ²    vs    Var(X_Poisson) = μ
            #
            # This makes the NB test more conservative (harder to call
            # enrichment) but more calibrated: when the null is true, the
            # fraction of false positives matches the stated α threshold.
            #
            # scipy.stats.nbinom uses a different parameterisation from the
            # (μ, α) form we reason about.  The conversion is:
            #
            #   n  = 1 / α          ("number of successes" in scipy's NB)
            #   p  = n / (n + μ)    ("probability of success" in scipy's NB)
            #
            # Deriving this:
            #   scipy's NB has mean = n(1-p)/p and variance = n(1-p)/p²
            #   Setting mean = μ:   n(1-p)/p = μ   →   n = μ·p/(1-p)
            #   Setting var  = μ + α·μ²:
            #     n(1-p)/p² = μ + α·μ²
            #     μ/p       = μ + α·μ²     (substituting mean formula)
            #     1/p       = 1 + α·μ
            #     p         = 1/(1 + α·μ)  = n/(n + μ)  where n = 1/α
            #
            # So: n_nb = 1/α,  p_nb = n_nb / (n_nb + expected)

            n_nb = 1.0 / nb_alpha
            p_nb = n_nb / (n_nb + expected)

            # One-sided p-value: P(X >= cnt | NB(n_nb, p_nb))
            #
            # nbinom.cdf(k, n, p) = P(X <= k)
            # Therefore P(X >= cnt) = 1 - P(X <= cnt-1) = 1 - nbinom.cdf(cnt-1, n, p)
            #
            # Using cnt-1 (not cnt) is essential for discrete distributions:
            # it includes the observed value itself in the tail being tested.
            # If we used cnt, we would compute P(X > cnt) which excludes the
            # observed value — a logical error for a test of "at least as extreme".
            pval_nb = 1 - stats.nbinom.cdf(cnt - 1, n=n_nb, p=p_nb)

        else:
            # If anchor data or region size is missing, set to NaN so these
            # rows are excluded from FDR correction rather than distorting it.
            expected        = np.nan
            fold_enrichment = np.nan
            pval_nb    = np.nan

        rows.append({
            "sample_id":       sid,
            "cell_type":       ct,
            "raw_count":       cnt,
            "norm_count":      norm_cnt,
            "zscore":          z,
            "expected_count":  expected,
            "fold_enrichment": fold_enrichment,
            "pval_nb":         pval_nb,   # column renamed to reflect NB test
        })

# ── Build DataFrame and apply BH FDR correction within each sample ───────────
# We correct across all cell types tested within a single sample (up to 16
# tests per sample).  Correcting within-sample rather than globally keeps the
# FDR meaningful: we are controlling the false discovery rate among cell types
# called enriched in that individual sample.
scores_df = pd.DataFrame(rows).merge(meta[["sample_id", "group"]], on="sample_id")


def _add_fdr(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Apply BH correction to pval_nb within one sample's cell types."""
    pvals = sample_df["pval_nb"].values
    # Only correct non-NaN values; NaN rows keep padj = NaN
    valid_mask = ~np.isnan(pvals)
    padj = np.full(len(pvals), np.nan)
    if valid_mask.sum() > 1:
        # Multiple tests: apply BH
        padj[valid_mask] = _bh_correct(pvals[valid_mask])
    elif valid_mask.sum() == 1:
        # Single valid test: padj equals pval (no correction needed)
        padj[valid_mask] = pvals[valid_mask]
    sample_df = sample_df.copy()
    sample_df["padj_nb"] = padj
    # enriched: True if the BH-adjusted p-value is below 0.05.
    # This is the binary per-sample per-cell-type enrichment call.
    sample_df["enriched"] = sample_df["padj_nb"] < 0.05
    return sample_df


# groupby sample_id applies _add_fdr independently to each sample's rows
scores_df = scores_df.groupby("sample_id", group_keys=False).apply(_add_fdr)

# ─── Step 6: Composite indices (CDI, NII, IAI) ───────────────────────────────
pivot = scores_df.pivot(index="sample_id", columns="cell_type", values="zscore")

composite = pd.DataFrame({"sample_id": list(pivot.index)})
for name, types in [("CDI", CDI_TYPES), ("NII", NII_TYPES), ("IAI", IAI_TYPES)]:
    cols = [c for c in types if c in pivot.columns]
    composite[name] = pivot[cols].mean(axis=1).values if cols else np.nan
composite = composite.merge(meta[["sample_id", "group"]], on="sample_id")

# ─── Step 7: Save TSVs ───────────────────────────────────────────────────────
Path(snakemake.output.scores).parent.mkdir(parents=True, exist_ok=True)
scores_df.to_csv(snakemake.output.scores, sep="\t", index=False)
composite.to_csv(snakemake.output.composite, sep="\t", index=False)

# ─── Step 8: Heatmap ─────────────────────────────────────────────────────────
group_order = ["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
group_colors = {
    "Ctrl": "#4878CF", "NEW": "#6ACC65",
    "MS-Rituximab-Stable": "#D65F5F",
    "MS-Rituximab-Progressive": "#B47CC7",
}

_mi = meta.set_index("sample_id")
sample_order = (
    _mi.loc[_mi["group"].isin(group_order), "group"]
    .sort_values()
    .index.tolist()
)
sample_order = [s for s in sample_order if s in pivot.index]

fig, ax = plt.subplots(
    figsize=(max(8, len(sample_order) * 0.18), max(6, len(pivot.columns) * 0.4))
)
sns.heatmap(
    pivot.loc[sample_order].T,
    cmap="RdBu_r", center=0, vmin=-3, vmax=3,
    ax=ax, yticklabels=True, xticklabels=False,
    cbar_kws={"label": "Poisson z-score"},
)

grp_colors = [
    group_colors.get(_mi.loc[s, "group"], "grey")
    for s in sample_order if s in _mi.index
]
for i, c in enumerate(grp_colors):
    ax.add_patch(plt.Rectangle(
        (i, len(pivot.columns)), 1, 0.5,
        color=c, clip_on=False, transform=ax.transData,
    ))

ax.set_title("Cell-type deconvolution (Poisson z-scores)")
ax.set_xlabel("Samples")
plt.tight_layout()
fig.savefig(snakemake.output.heatmap)
plt.close(fig)

# ─── Step 9: Violin plots per cell type with between-group statistics ─────────
group_palette = {
    "Ctrl": "#4878CF", "NEW": "#6ACC65",
    "MS-Rituximab-Stable": "#D65F5F",
    "MS-Rituximab-Progressive": "#B47CC7",
}

def dunn_posthoc(df, group_col, value_col, groups):
    """Dunn's test approximation: Mann-Whitney U with Bonferroni correction."""
    pairs = list(combinations(groups, 2))
    results = {}
    for g1, g2 in pairs:
        a = df.loc[df[group_col] == g1, value_col].dropna()
        b = df.loc[df[group_col] == g2, value_col].dropna()
        if len(a) < 2 or len(b) < 2:
            results[(g1, g2)] = 1.0
            continue
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        results[(g1, g2)] = p
    # Bonferroni correction
    n = len(pairs)
    return {k: min(v * n, 1.0) for k, v in results.items()}


def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


present_groups = [g for g in group_order if g in scores_df["group"].unique()]
cell_types_list = list(CELL_TYPES.keys())

# ─── Step 10: Collect all figures and write single PDF ───────────────────────
# Two new figures appended to the violin PDF.
#
# Figure A — Fold-enrichment heatmap
#   Same layout as the z-score heatmap (Step 8) but the colour encodes
#   fold_enrichment = observed / expected, where expected comes from the
#   within-sample anchor density.  This is group-agnostic: a value of 2
#   means "this sample has twice as many midpoints at this cell type's
#   signature than the constitutive baseline predicts", regardless of what
#   the controls look like.  Asterisks mark cells where padj_nb < 0.05.
#
# Figure B — Enrichment rate bar chart
#   For each cell type, one grouped bar per disease group showing the
#   fraction of samples in that group called enriched (padj < 0.05).
#   Fisher's exact test is run between every pair of groups; significant
#   pairwise comparisons (p < 0.05) are annotated with brackets.

from scipy.stats import fisher_exact as _fisher_exact

# ── Figure A: fold-enrichment heatmap ────────────────────────────────────────

# Pivot fold_enrichment into a samples × cell-types matrix, same ordering as
# the z-score pivot used above.
fe_pivot = scores_df.pivot(
    index="sample_id", columns="cell_type", values="fold_enrichment"
)
# Pivot enriched flag into a parallel boolean matrix for asterisk annotation.
enriched_pivot = scores_df.pivot(
    index="sample_id", columns="cell_type", values="enriched"
).fillna(False)

# Align both pivots to the same sample and cell-type order.
fe_sample_order = [s for s in sample_order if s in fe_pivot.index]
fe_ct_order     = [c for c in fe_pivot.columns]

# Collect all figures in a list; write to violin PDF in one PdfPages call
# at the very end so all pages (violin + composite + Poisson enrichment) land
# in a single file without needing to reopen/append.
_all_figs = []

# ── Re-generate violin figures (same logic as Step 9) so we can collect them --
for ct in cell_types_list:
    ct_df = scores_df[scores_df["cell_type"] == ct].copy()
    if ct_df.empty:
        continue
    fig, ax = plt.subplots(figsize=(7, 5))
    plot_df = ct_df[ct_df["group"].isin(present_groups)].copy()
    plot_df["group"] = pd.Categorical(plot_df["group"], categories=present_groups)
    plot_df = plot_df.sort_values("group")
    group_data = [plot_df.loc[plot_df["group"] == g, "zscore"].dropna().values
                  for g in present_groups]
    colors = [group_palette.get(g, "grey") for g in present_groups]
    parts = ax.violinplot(group_data, positions=range(len(present_groups)),
                          showmedians=True, showextrema=True)
    for pc, col in zip(parts["bodies"], colors):
        pc.set_facecolor(col); pc.set_alpha(0.6)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        parts[key].set_color("black")
    rng = np.random.default_rng(42)
    for i, (grp, col) in enumerate(zip(present_groups, colors)):
        vals = plot_df.loc[plot_df["group"] == grp, "zscore"].dropna()
        jx   = rng.uniform(-0.08, 0.08, size=len(vals)) + i
        ax.scatter(jx, vals, color=col, s=20, zorder=3, alpha=0.8,
                   edgecolors="white", lw=0.4)
    ax.set_xticks(range(len(present_groups)))
    ax.set_xticklabels([g.replace("MS-Rituximab-", "RTX-") for g in present_groups],
                       rotation=20, ha="right", fontsize=9)
    ax.axhline(0, ls="--", lw=0.8, color="grey")
    ax.set_ylabel("Poisson z-score vs Ctrl")
    ax.set_title(ct)
    valid = [d for d in group_data if len(d) >= 2]
    if len(valid) >= 2:
        kw_stat, kw_p = stats.kruskal(*valid)
        ax.text(0.98, 0.98, f"KW p={kw_p:.3g}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8)
        posthoc  = dunn_posthoc(plot_df, "group", "zscore", present_groups)
        sig_pairs = {k: v for k, v in posthoc.items() if v < 0.05}
        y_top = plot_df["zscore"].max() if not plot_df.empty else 3
        step  = (y_top - ax.get_ylim()[0]) * 0.12
        for lvl, ((g1, g2), p) in enumerate(sig_pairs.items()):
            x1, x2 = present_groups.index(g1), present_groups.index(g2)
            y = y_top + step * (lvl + 1)
            ax.plot([x1, x1, x2, x2], [y - step*0.2, y, y, y - step*0.2],
                    lw=1, color="black")
            ax.text((x1+x2)/2, y, sig_label(p), ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    _all_figs.append(fig)

for idx_col in ["CDI", "NII", "IAI"]:
    if idx_col not in composite.columns:
        continue
    fig, ax = plt.subplots(figsize=(7, 5))
    comp_sub = composite[composite["group"].isin(present_groups)].copy()
    comp_sub["group"] = pd.Categorical(comp_sub["group"], categories=present_groups)
    group_data = [comp_sub.loc[comp_sub["group"] == g, idx_col].dropna().values
                  for g in present_groups]
    colors = [group_palette.get(g, "grey") for g in present_groups]
    parts = ax.violinplot(group_data, positions=range(len(present_groups)),
                          showmedians=True, showextrema=True)
    for pc, col in zip(parts["bodies"], colors):
        pc.set_facecolor(col); pc.set_alpha(0.6)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        parts[key].set_color("black")
    rng2 = np.random.default_rng(0)
    for i, (grp, col) in enumerate(zip(present_groups, colors)):
        vals = comp_sub.loc[comp_sub["group"] == grp, idx_col].dropna()
        jx   = rng2.uniform(-0.08, 0.08, size=len(vals)) + i
        ax.scatter(jx, vals, color=col, s=22, zorder=3, alpha=0.8,
                   edgecolors="white", lw=0.4)
    ax.set_xticks(range(len(present_groups)))
    ax.set_xticklabels([g.replace("MS-Rituximab-", "RTX-") for g in present_groups],
                       rotation=20, ha="right", fontsize=9)
    ax.axhline(0, ls="--", lw=0.8, color="grey")
    ax.set_ylabel("Composite index score")
    ax.set_title(f"{idx_col} — composite index")
    valid = [d for d in group_data if len(d) >= 2]
    if len(valid) >= 2:
        _, kw_p = stats.kruskal(*valid)
        ax.text(0.98, 0.98, f"KW p={kw_p:.3g}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8)
    plt.tight_layout()
    _all_figs.append(fig)

# ── Figure A: fold-enrichment heatmap ────────────────────────────────────────
if len(fe_sample_order) > 0 and len(fe_ct_order) > 0:
    fig, ax = plt.subplots(
        figsize=(max(8, len(fe_sample_order) * 0.18),
                 max(6, len(fe_ct_order) * 0.4))
    )
    # log2(fold_enrichment) centres the colour scale at 0 (fold=1 → log2=0),
    # making enrichment (>1) red and depletion (<1) blue, matching the z-score
    # heatmap colour logic.  Clip to avoid log2(0) = -inf.
    log2_fe = np.log2(
        fe_pivot.loc[fe_sample_order, fe_ct_order]
        .clip(lower=0.01)
        .astype(float)
    )
    sns.heatmap(
        log2_fe.T,
        cmap="RdBu_r", center=0, vmin=-3, vmax=3,
        ax=ax, yticklabels=True, xticklabels=False,
        cbar_kws={"label": "log2(fold enrichment over anchor baseline)"},
    )
    # Overlay an asterisk on cells where padj_nb < 0.05.
    # enriched_pivot rows = samples, columns = cell types; we iterate over
    # the transposed layout (rows = cell types, columns = samples).
    enr_mat = enriched_pivot.loc[fe_sample_order, fe_ct_order].T.values
    for row_i, ct_name in enumerate(fe_ct_order):
        for col_j, sid in enumerate(fe_sample_order):
            if enr_mat[row_i, col_j]:
                # Place asterisk at the centre of the cell
                ax.text(col_j + 0.5, row_i + 0.5, "*",
                        ha="center", va="center",
                        fontsize=6, color="black", fontweight="bold")
    # Group colour bar below x-axis (same as z-score heatmap)
    fe_grp_colors = [
        group_colors.get(_mi.loc[s, "group"], "grey")
        for s in fe_sample_order if s in _mi.index
    ]
    for i, c in enumerate(fe_grp_colors):
        ax.add_patch(plt.Rectangle(
            (i, len(fe_ct_order)), 1, 0.5,
            color=c, clip_on=False, transform=ax.transData,
        ))
    ax.set_title("Cell-type deconvolution — log2 fold enrichment over anchor baseline\n"
                 "(* = padj_nb < 0.05 within sample)")
    ax.set_xlabel("Samples")
    plt.tight_layout()
    _all_figs.append(fig)

# ── Figure B: enrichment rate bar chart with Fisher's exact test ─────────────
# For each cell type: one bar per group = fraction of samples called enriched.
# Fisher's exact test on a 2×2 contingency table (enriched vs not, group1 vs
# group2) tests whether the enrichment rate differs between groups.
#
# Contingency table layout:
#              enriched    not enriched
#   group 1  [  n1_yes       n1_no   ]
#   group 2  [  n2_yes       n2_no   ]
#
# fisher_exact returns (odds_ratio, p_value); we use p_value only.

for ct in cell_types_list:
    ct_data = scores_df[scores_df["cell_type"] == ct].copy()
    if ct_data.empty or "enriched" not in ct_data.columns:
        continue

    # Compute enrichment rate (fraction True) and sample size per group
    rates   = {}
    n_per   = {}
    for grp in present_groups:
        grp_vals = ct_data.loc[ct_data["group"] == grp, "enriched"].dropna()
        rates[grp] = float(grp_vals.mean()) if len(grp_vals) > 0 else 0.0
        n_per[grp] = len(grp_vals)

    fig, ax = plt.subplots(figsize=(7, 5))
    x      = np.arange(len(present_groups))
    colors = [group_palette.get(g, "grey") for g in present_groups]

    # Draw bars: height = fraction enriched, label shows n above bar
    bars = ax.bar(x, [rates[g] for g in present_groups],
                  color=colors, alpha=0.8, edgecolor="black", linewidth=0.6)
    for bar, grp in zip(bars, present_groups):
        n = n_per[grp]
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.02,
                f"n={n}", ha="center", va="bottom", fontsize=8)

    ax.set_ylim(0, 1.25)
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace("MS-Rituximab-", "RTX-") for g in present_groups],
                       rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Fraction of samples enriched (padj < 0.05)")
    ax.set_title(f"{ct} — enrichment rate by group")
    ax.axhline(0.05, ls=":", lw=0.8, color="grey", label="5% reference")

    # Fisher's exact test between every pair of groups;
    # annotate only significant pairs (p < 0.05) with brackets.
    pair_results = []
    for g1, g2 in combinations(present_groups, 2):
        d1 = ct_data.loc[ct_data["group"] == g1, "enriched"].dropna()
        d2 = ct_data.loc[ct_data["group"] == g2, "enriched"].dropna()
        if len(d1) < 2 or len(d2) < 2:
            continue
        # Build 2×2 table: [[enriched_g1, not_g1], [enriched_g2, not_g2]]
        n1_yes = int(d1.sum());  n1_no = len(d1) - n1_yes
        n2_yes = int(d2.sum());  n2_no = len(d2) - n2_yes
        _, p_fisher = _fisher_exact([[n1_yes, n1_no], [n2_yes, n2_no]])
        if p_fisher < 0.05:
            pair_results.append((g1, g2, p_fisher))

    # Draw significance brackets above the bars
    y_base = max(rates.values()) + 0.12 if rates else 0.8
    for lvl, (g1, g2, p_fisher) in enumerate(pair_results):
        x1 = present_groups.index(g1)
        x2 = present_groups.index(g2)
        y  = y_base + lvl * 0.12
        ax.plot([x1, x1, x2, x2], [y - 0.02, y, y, y - 0.02],
                lw=1, color="black")
        ax.text((x1 + x2) / 2, y + 0.01, sig_label(p_fisher),
                ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    _all_figs.append(fig)

# ── Save all figures into the violin PDF in one pass ─────────────────────────
with pdf_backend.PdfPages(snakemake.output.violin) as _pp:
    for _fig in _all_figs:
        _pp.savefig(_fig)
        plt.close(_fig)

print("[Module 4] Done.")
