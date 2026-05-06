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

import gzip
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
import statsmodels.api as sm


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


def _count_per_region_sample(args):
    """
    Per-region midpoint counting for the within-sample enrichment analysis.
    Unlike _count_one_sample (which sums by cell-type label), this function
    returns one integer count per row of combined_reg_bed.  col4 of that BED
    must hold a unique region_id per row.
    Returns: (sample_id, {region_id: int_count})
    """
    sid, frags_bed, combined_reg_bed = args
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
            f"bedtools intersect -a {combined_reg_bed} -b {mid_tmp} -c",
            shell=True, capture_output=True, text=True,
        )
    finally:
        if os.path.exists(mid_tmp):
            os.remove(mid_tmp)
    per_region: dict[str, int] = {}
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 5:
            per_region[parts[3]] = int(parts[-1])
    return sid, per_region


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

# === within-sample enrichment (group-agnostic) — Sadeh local background ===
#
# Background = 5 kb tiles from the complement of (RRE ∪ DAC) on autosomes.
# Every tile is guaranteed non-overlapping with signal regions and the non-RRE
# genome is covered contiguously (bedtools complement + makewindows -s = -w).
# Per-sample expected density fitted at three levels (genome → chr → 2 Mb),
# mirroring Sadeh 2021 Background.R buildBackground().
# Each sample is scored independently — no cohort statistic is needed.

# ── WS-A: constants ───────────────────────────────────────────────────────────
_WS_BG_WIN   = 5_000        # bp — background tile size
_WS_REGION   = 2_000_000    # bp — neighbourhood for local rate fitting
_WS_BG_TRIM  = 0.95         # Sadeh fitNoise percentile trim
_WS_MIN_BG   = 20           # min windows to fit a regional rate
_WS_MIN_SIG  = 50           # min signature windows for NB GLM
_WS_N_PERM   = 200          # permuted-null iterations
_WS_PERM_THR = 0.2          # |median_perm_log2FE| flag threshold
_WS_EPS      = 1.0          # log2FE pseudocount

# ── WS-B: build background windows (complement of RRE+DAC, tiled at 5 kb) ────
print("[Module 4] Within-sample enrichment: building background windows ...")

_ws_chrom_sizes = snakemake.params.chrom_sizes
_ws_dac         = snakemake.input.dac_regions
_ws_rre_bed     = snakemake.input.rre_universe

# Autosomal chrom.sizes sorted lexicographically (same order as `sort -k1,1`
# used on the merged BED — bedtools complement requires matching sort orders).
_ws_t_auto_unsorted = tempfile.NamedTemporaryFile(
    suffix=".txt", delete=False, mode="w")
with open(_ws_chrom_sizes) as _fh_cs:
    for _line_cs in _fh_cs:
        _p_cs = _line_cs.split()
        if _p_cs and _p_cs[0].startswith("chr") and _p_cs[0][3:].isdigit():
            _ws_t_auto_unsorted.write(_line_cs)
_ws_t_auto_unsorted.close()
_ws_auto_sizes = tempfile.NamedTemporaryFile(
    suffix=".txt", delete=False).name
subprocess.run(
    f"sort -k1,1 {_ws_t_auto_unsorted.name} > {_ws_auto_sizes}",
    shell=True, check=True)
os.remove(_ws_t_auto_unsorted.name)

_ws_t_merged = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
_ws_t_compl  = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
_ws_t_tiled  = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
_ws_t_bg     = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name

# a. Merge RRE + DAC on autosomes; clip to chrom boundaries so complement
#    never produces negative-width intervals (RRE/DAC may extend beyond
#    the chromosome end listed in the genome file).
subprocess.run(
    f"cat {_ws_rre_bed} {_ws_dac}"
    f" | grep -P '^chr[0-9]+\\t'"
    f" | cut -f1-3"
    f" | sort -k1,1 -k2,2n"
    f" | bedtools merge"
    f" | bedtools intersect -a - -b <(awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,0,$2}}'"
    f" {_ws_auto_sizes}) -u"
    f" > {_ws_t_merged}",
    shell=True, check=True, executable="/bin/bash")

# b. Complement on autosomes; filter degenerate intervals (start >= end)
#    that can arise when a merged region touches the chromosome boundary.
subprocess.run(
    f"bedtools complement -i {_ws_t_merged} -g {_ws_auto_sizes}"
    f" | awk '$3 > $2'"
    f" > {_ws_t_compl}",
    shell=True, check=True)

# c. Tile with no overlap (step = window → non-overlapping, contiguous)
subprocess.run(
    f"bedtools makewindows -b {_ws_t_compl}"
    f" -w {_WS_BG_WIN} -s {_WS_BG_WIN}"
    f" > {_ws_t_tiled}",
    shell=True, check=True)

# d. Drop partial end-of-interval tiles
subprocess.run(
    f"awk '$3-$2=={_WS_BG_WIN}' {_ws_t_tiled} > {_ws_t_bg}",
    shell=True, check=True)

for _f_tmp in [_ws_t_merged, _ws_t_compl, _ws_t_tiled, _ws_auto_sizes]:
    if os.path.exists(_f_tmp):
        os.remove(_f_tmp)

_ws_bg_df = pd.read_csv(_ws_t_bg, sep="\t", header=None,
                        names=["chr", "start", "end"])
os.remove(_ws_t_bg)

_ws_bg_df["region_id"]  = ("bg|" + _ws_bg_df["chr"] + ":"
                            + _ws_bg_df["start"].astype(str) + "-"
                            + _ws_bg_df["end"].astype(str))
_ws_bg_df["length_kb"]  = (_ws_bg_df["end"] - _ws_bg_df["start"]) / 1000.0
_ws_bg_df["region_2mb"] = (_ws_bg_df["chr"] + ":"
                            + (_ws_bg_df["start"] // _WS_REGION).astype(str))

if len(_ws_bg_df) < 2000:
    raise RuntimeError(
        f"[Module 4] ABORT: only {len(_ws_bg_df)} background windows after filtering. "
        "Check RRE universe, DAC, and chrom sizes.")
print(f"[Module 4]   {len(_ws_bg_df):,} background windows.")


# ── WS-C: per-cell-type signal region metadata ────────────────────────────────
# ct_regions[ct] has integer column keys: 0=chr, 1=start, 2=end
# Signal region_id: "sig|{ct}|{chr}:{start}-{end}"

_ws_ct_sig: dict[str, pd.DataFrame] = {}
for _ct_ws in CELL_TYPES:
    _df_ct = ct_regions.get(_ct_ws, pd.DataFrame()).copy()
    if _df_ct.empty:
        _ws_ct_sig[_ct_ws] = _df_ct
        continue
    _df_ct = _df_ct.reset_index(drop=True)
    _df_ct["region_id"] = ("sig|" + _ct_ws + "|"
                            + _df_ct[0].astype(str) + ":"
                            + _df_ct[1].astype(int).astype(str) + "-"
                            + _df_ct[2].astype(int).astype(str))
    _df_ct["length_kb"]  = (_df_ct[2] - _df_ct[1]).astype(float) / 1000.0
    _df_ct["chr"]        = _df_ct[0].astype(str)
    _df_ct["start"]      = _df_ct[1].astype(int)
    _df_ct["region_2mb"] = (_df_ct["chr"] + ":"
                             + (_df_ct["start"] // _WS_REGION).astype(str))
    _ws_ct_sig[_ct_ws] = _df_ct

# ── WS-D: build combined counting BED and save gzipped copy ──────────────────
print("[Module 4] Building combined counting BED ...")

_ws_sig_rows: list[str] = []
for _ct_ws, _df_ct in _ws_ct_sig.items():
    for _, _row_ct in _df_ct.iterrows():
        _ws_sig_rows.append(
            f"{_row_ct['chr']}\t{_row_ct['start']}\t"
            f"{int(_row_ct[2])}\t{_row_ct['region_id']}")

_ws_bg_rows: list[str] = [
    f"{_r.chr}\t{_r.start}\t{_r.end}\t{_r.region_id}"
    for _, _r in _ws_bg_df.iterrows()
]

_ws_t_unsorted = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
_ws_t_unsorted.write("\n".join(_ws_sig_rows + _ws_bg_rows) + "\n")
_ws_t_unsorted.close()

_ws_combined = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
subprocess.run(
    f"sort -k1,1 -k2,2n {_ws_t_unsorted.name} > {_ws_combined}",
    shell=True, check=True)
os.remove(_ws_t_unsorted.name)

with open(_ws_combined, "rb") as _fh_comb, \
     gzip.open(snakemake.output.within_sample_windows_bed, "wb") as _fh_gz:
    _fh_gz.write(_fh_comb.read())
print(f"[Module 4]   {len(_ws_sig_rows):,} signal + {len(_ws_bg_rows):,} bg rows "
      f"→ {snakemake.output.within_sample_windows_bed}")

# ── WS-E: parallel fragment counting ─────────────────────────────────────────
print(f"[Module 4] Per-region counting ({len(samples)} samples, {n_cpus} workers) ...")
_ws_count_args = [(sid, f"{frags_dir}/{sid}.bed", _ws_combined) for sid in samples]
_ws_sample_counts: dict[str, dict[str, int]] = {}
with ProcessPoolExecutor(max_workers=n_cpus) as _ws_pool:
    for _sid_ws, _cnts_ws in _ws_pool.map(_count_per_region_sample, _ws_count_args):
        _ws_sample_counts[_sid_ws] = _cnts_ws
        print(f"  done: {_sid_ws}")
os.remove(_ws_combined)
print("[Module 4] Counting complete.")

# ── WS-F: per-sample Sadeh background model ───────────────────────────────────
_ws_bg_rids   = _ws_bg_df["region_id"].values
_ws_bg_len_kb = _ws_bg_df["length_kb"].values
_ws_bg_chr    = _ws_bg_df["chr"].values
_ws_bg_start  = _ws_bg_df["start"].values
_ws_bg_reg2mb = _ws_bg_df["region_2mb"].values


def _ws_fit_noise(densities: np.ndarray) -> float:
    """95th-pct-trimmed mean density — mirrors Sadeh Background.R fitNoise()."""
    x = densities[np.isfinite(densities) & (densities >= 0)]
    if len(x) < _WS_MIN_BG:
        return float("nan")
    x = x[x <= np.quantile(x, _WS_BG_TRIM)]
    return float(np.mean(x)) if len(x) > 0 else float("nan")


def _ws_build_bg_model(sid: str) -> dict:
    cnts    = _ws_sample_counts[sid]
    counts  = np.array([cnts.get(r, 0) for r in _ws_bg_rids], dtype=float)
    density = counts / _ws_bg_len_kb
    rate_g  = _ws_fit_noise(density)
    rate_chr: dict[str, float] = {}
    for _ch in np.unique(_ws_bg_chr):
        _r = _ws_fit_noise(density[_ws_bg_chr == _ch])
        rate_chr[_ch] = _r if np.isfinite(_r) else rate_g
    rate_reg: dict[str, float] = {}
    for _rk in np.unique(_ws_bg_reg2mb):
        _r  = _ws_fit_noise(density[_ws_bg_reg2mb == _rk])
        _ch = _rk.rsplit(":", 1)[0]
        _fb = rate_chr.get(_ch, rate_g)
        rate_reg[_rk] = _r if np.isfinite(_r) else _fb
    return {"genome": rate_g, "chr": rate_chr, "region": rate_reg}


def _ws_expected_density(chrom: str, start: int, model: dict) -> float:
    _rk = f"{chrom}:{start // _WS_REGION}"
    _r  = model["region"].get(_rk)
    if _r is None or not np.isfinite(_r):
        _r = model["chr"].get(chrom)
    if _r is None or not np.isfinite(_r):
        _r = model["genome"]
    return max(float(_r), 1e-6) if (_r is not None and np.isfinite(_r)) else 1e-6


print("[Module 4] Fitting per-sample background models ...")
_ws_bg_models = {sid: _ws_build_bg_model(sid)
                 for sid in samples if sid in _ws_sample_counts}

_ws_bg_chr_idx: dict[str, np.ndarray] = {
    _ch: np.where(_ws_bg_chr == _ch)[0]
    for _ch in np.unique(_ws_bg_chr)
}

# ── WS-G: NB GLM helper ───────────────────────────────────────────────────────
def _ws_nb_glm(sig_cnt, bg_cnt, sig_len_kb, bg_len_kb):
    counts  = np.r_[sig_cnt, bg_cnt].astype(float)
    lengths = np.maximum(np.r_[sig_len_kb, bg_len_kb], 1e-3)
    is_sig  = np.r_[np.ones(len(sig_cnt)), np.zeros(len(bg_cnt))]
    X = sm.add_constant(is_sig, has_constant="add")
    try:
        fit = sm.GLM(counts, X, family=sm.families.NegativeBinomial(),
                     offset=np.log(lengths)).fit(disp=False, maxiter=100)
        return float(fit.params[1]), float(fit.pvalues[1])
    except Exception:
        return float("nan"), float("nan")

# ── WS-H: main enrichment loop ────────────────────────────────────────────────
print("[Module 4] Computing within-sample enrichment scores ...")
_ws_rng   = np.random.default_rng(42)
_ws_rows: list[dict] = []

for sid in samples:
    if sid not in _ws_sample_counts:
        continue
    _cnts_s  = _ws_sample_counts[sid]
    _model_s = _ws_bg_models[sid]

    for _ct_ws, _df_ct in _ws_ct_sig.items():
        if len(_df_ct) == 0:
            continue

        _sig_cnt = np.array([_cnts_s.get(_rid, 0)
                              for _rid in _df_ct["region_id"]], dtype=float)
        _sig_len = _df_ct["length_kb"].values

        _exp_dens = np.array(
            [_ws_expected_density(_r["chr"], _r["start"], _model_s)
             for _, _r in _df_ct.iterrows()], dtype=float)

        _sig_dens = _sig_cnt / np.maximum(_sig_len, 1e-3)
        _med_sig  = float(np.median(_sig_dens))
        _med_exp  = float(np.median(_exp_dens))

        _log2fe = float(np.log2((_med_sig + _WS_EPS) / (_med_exp + _WS_EPS)))
        _excess = max(0.0, _med_sig - _med_exp)

        _low_n   = len(_df_ct) < _WS_MIN_SIG
        _beta_nb = _p_nb = float("nan")
        if not _low_n:
            _bg_idx = np.concatenate(
                [_ws_bg_chr_idx.get(_ch, np.array([], int))
                 for _ch in _df_ct["chr"].unique()])
            if len(_bg_idx) > 0:
                _chosen = _ws_rng.choice(
                    _bg_idx,
                    size=min(len(_bg_idx), len(_df_ct) * 10),
                    replace=False)
                _beta_nb, _p_nb = _ws_nb_glm(
                    _sig_cnt,
                    np.array([_cnts_s.get(_ws_bg_rids[i], 0)
                               for i in _chosen], float),
                    _sig_len,
                    _ws_bg_len_kb[_chosen])

        _ws_rows.append({
            "sample": sid, "cell_type": _ct_ws,
            "log2FE": _log2fe, "beta_nb": _beta_nb, "p_nb": _p_nb,
            "q_nb": float("nan"), "excess": _excess,
            "contribution": float("nan"), "enriched": False,
            "low_n": _low_n, "n_sig_windows": len(_df_ct),
        })

ws_df_out = pd.DataFrame(_ws_rows)


def _ws_bh_fdr(grp: pd.DataFrame) -> pd.DataFrame:
    from statsmodels.stats.multitest import multipletests as _mtest
    pvals = grp["p_nb"].values
    valid = np.isfinite(pvals)
    q     = np.full(len(pvals), float("nan"))
    if valid.sum() > 1:
        _, q_vals, _, _ = _mtest(pvals[valid], method="fdr_bh")
        q[valid] = q_vals
    grp = grp.copy()
    grp["q_nb"]     = q
    grp["enriched"] = (q < 0.05) & (grp["beta_nb"].values > 0)
    return grp


ws_df_out = ws_df_out.groupby("sample", group_keys=False).apply(_ws_bh_fdr)


def _ws_contrib(grp: pd.DataFrame) -> pd.DataFrame:
    exc   = grp["excess"].clip(lower=0).fillna(0)
    total = exc.sum()
    grp   = grp.copy()
    grp["contribution"] = (exc / total) if total > 0 else 0.0
    return grp


ws_df_out = ws_df_out.groupby("sample", group_keys=False).apply(_ws_contrib)

# ── WS-I: sanity checks ───────────────────────────────────────────────────────
print("[Module 4] Running sanity checks ...")
_ws_perm_size = max(50, len(_ws_bg_rids) // 20)
_ws_qc_rows: list[dict] = []

for sid in samples:
    if sid not in _ws_sample_counts:
        continue
    _cnts_s  = _ws_sample_counts[sid]
    _model_s = _ws_bg_models[sid]

    # Permuted null: random background windows treated as fake signal
    _perm_fes: list[float] = []
    for _ in range(_WS_N_PERM):
        _idx_p  = _ws_rng.integers(0, len(_ws_bg_rids), size=_ws_perm_size)
        _p_cnt  = np.array([_cnts_s.get(_ws_bg_rids[i], 0)
                             for i in _idx_p], float)
        _p_dens = _p_cnt / _ws_bg_len_kb[_idx_p]
        _p_exp  = np.array(
            [_ws_expected_density(_ws_bg_chr[i], int(_ws_bg_start[i]), _model_s)
             for i in _idx_p], float)
        _perm_fes.append(float(np.log2(
            (np.median(_p_dens) + _WS_EPS) / (np.median(_p_exp) + _WS_EPS))))
    _med_perm = float(np.median(_perm_fes))

    # Anti-signal: separate random set of background windows
    _anti_idx    = _ws_rng.integers(0, len(_ws_bg_rids),
                                    size=min(500, len(_ws_bg_rids)))
    _a_dens      = (np.array([_cnts_s.get(_ws_bg_rids[i], 0)
                               for i in _anti_idx], float)
                    / _ws_bg_len_kb[_anti_idx])
    _a_exp       = np.array(
        [_ws_expected_density(_ws_bg_chr[i], int(_ws_bg_start[i]), _model_s)
         for i in _anti_idx], float)
    _anti_log2fe = float(np.log2(
        (np.median(_a_dens) + _WS_EPS) / (np.median(_a_exp) + _WS_EPS)))

    _n_enr = int(((ws_df_out["sample"] == sid) & ws_df_out["enriched"]).sum())
    _ws_qc_rows.append({
        "sample": sid,
        "median_perm_log2FE":    _med_perm,
        "neg_ctrl_log2FE":       _anti_log2fe,
        "n_signatures_enriched": _n_enr,
        "perm_flag":             abs(_med_perm) > _WS_PERM_THR,
    })

qc_ws_df = pd.DataFrame(_ws_qc_rows)
_n_perm_flagged = int(qc_ws_df["perm_flag"].sum())
if _n_perm_flagged > 0:
    print(f"[Module 4] WARNING: {_n_perm_flagged}/{len(qc_ws_df)} samples "
          f"flagged by permuted null (|median_perm_log2FE| > {_WS_PERM_THR}).")
if _n_perm_flagged == len(samples) and len(samples) > 0:
    raise RuntimeError(
        "[Module 4] ABORT: permuted-signature null fails for every sample. "
        "Check RRE/DAC filters and chrom sizes.")

# ── WS-J: write TSV outputs ───────────────────────────────────────────────────
_ws_out_dir = Path(snakemake.output.within_sample_scores).parent
_ws_out_dir.mkdir(parents=True, exist_ok=True)

ws_df_out[[
    "sample", "cell_type", "log2FE", "beta_nb", "p_nb", "q_nb",
    "excess", "contribution", "enriched", "low_n", "n_sig_windows",
]].to_csv(snakemake.output.within_sample_scores, sep="\t", index=False)

qc_ws_df.to_csv(snakemake.output.within_sample_qc, sep="\t", index=False)
print(f"[Module 4] Scores → {snakemake.output.within_sample_scores}")

# ── WS-K: figures ─────────────────────────────────────────────────────────────
_ws_ct_order = list(CELL_TYPES.keys())

# Fig log2FE: group-agnostic boxplot (page 1) + group-coloured scatter (page 2)
_ws_box_data = [
    ws_df_out.loc[ws_df_out["cell_type"] == ct, "log2FE"].dropna().values
    for ct in _ws_ct_order
]
_fig_log2fe_ws, _ax_log2fe = plt.subplots(
    figsize=(max(8, len(_ws_ct_order) * 0.7), 5))
_ax_log2fe.boxplot(
    _ws_box_data, labels=_ws_ct_order,
    patch_artist=True, notch=False, sym=".",
    medianprops={"color": "black", "linewidth": 1.5},
    boxprops={"facecolor": "#AECDE8", "alpha": 0.8},
)
_ax_log2fe.axhline(0, ls="--", lw=0.8, color="grey")
_ax_log2fe.set_ylabel("log₂ Fold Enrichment (vs Sadeh local background)")
_ax_log2fe.set_title("Within-sample tissue enrichment (group-agnostic)")
_ax_log2fe.set_xticklabels(_ws_ct_order, rotation=45, ha="right", fontsize=8)
plt.tight_layout()

# Group-coloured variant (page 2 — downstream reference)
_ws_grp_df = ws_df_out.merge(
    meta[["sample_id", "group"]], left_on="sample", right_on="sample_id", how="left")
_n_ct_ws = len(_ws_ct_order)
_fig_grp_ws, _axes_grp = plt.subplots(
    1, _n_ct_ws, figsize=(max(12, _n_ct_ws * 1.1), 5), sharey=True)
if _n_ct_ws == 1:
    _axes_grp = [_axes_grp]
_ws_jitter_rng = np.random.default_rng(7)
for _ax_g, _ct_g in zip(_axes_grp, _ws_ct_order):
    _ct_sub = _ws_grp_df[_ws_grp_df["cell_type"] == _ct_g]
    for _grp_g, _col_g in group_palette.items():
        _vals_g = _ct_sub.loc[_ct_sub["group"] == _grp_g, "log2FE"].dropna()
        if len(_vals_g) > 0:
            _jx = _ws_jitter_rng.uniform(-0.15, 0.15, size=len(_vals_g))
            _ax_g.scatter(_jx, _vals_g, color=_col_g, s=12, alpha=0.7, label=_grp_g)
    _ax_g.axhline(0, ls="--", lw=0.6, color="grey")
    _ax_g.set_title(_ct_g, fontsize=7, rotation=45, ha="right")
    _ax_g.set_xticks([])
_handles_grp = [
    plt.Line2D([0], [0], marker="o", color=c, lw=0, label=g)
    for g, c in group_palette.items()
]
_fig_grp_ws.legend(handles=_handles_grp, loc="upper right", fontsize=7)
_fig_grp_ws.suptitle("Within-sample log2FE by group (Sadeh local background)")
plt.tight_layout()

with pdf_backend.PdfPages(snakemake.output.within_sample_log2fe) as _pp_ws:
    _pp_ws.savefig(_fig_log2fe_ws)  # page 1: group-agnostic
    _pp_ws.savefig(_fig_grp_ws)     # page 2: group-coloured variant
plt.close(_fig_log2fe_ws)
plt.close(_fig_grp_ws)

# Stacked bars: per-sample fractional contribution (Sadeh Fig 4c style)
_ws_ct_colors = sns.color_palette("tab20", n_colors=len(_ws_ct_order))
_ws_cmap      = dict(zip(_ws_ct_order, _ws_ct_colors))
_ws_contrib_piv = (ws_df_out
                   .pivot(index="sample", columns="cell_type", values="contribution")
                   .fillna(0).sort_index())

_fig_contrib_ws, _ax_contrib = plt.subplots(
    figsize=(max(10, len(_ws_contrib_piv) * 0.15), 5))
_ws_bottom = np.zeros(len(_ws_contrib_piv))
for _ct_cb in _ws_ct_order:
    if _ct_cb not in _ws_contrib_piv.columns:
        continue
    _vals_cb = _ws_contrib_piv[_ct_cb].values
    _ax_contrib.bar(range(len(_ws_contrib_piv)), _vals_cb, bottom=_ws_bottom,
                    color=_ws_cmap[_ct_cb], label=_ct_cb, width=1.0, edgecolor="none")
    _ws_bottom += _vals_cb
_ax_contrib.set_xlim(-0.5, len(_ws_contrib_piv) - 0.5)
_ax_contrib.set_ylim(0, 1)
_ax_contrib.set_ylabel("Fractional contribution")
_ax_contrib.set_title("Per-sample cell-type contributions (alphabetical order)")
_ax_contrib.set_xticks([])
_ax_contrib.set_xlabel("Samples (alphabetical, no group sorting)")
_ax_contrib.legend(loc="upper right", fontsize=6, ncol=2)
plt.tight_layout()
with pdf_backend.PdfPages(snakemake.output.within_sample_contribution) as _pp_cb:
    _pp_cb.savefig(_fig_contrib_ws)
plt.close(_fig_contrib_ws)

# QC scatter: permuted null vs anti-signal log2FE per sample
_fig_qc_ws, _ax_qc = plt.subplots(figsize=(6, 5))
_ax_qc.scatter(qc_ws_df["median_perm_log2FE"], qc_ws_df["neg_ctrl_log2FE"],
               s=30, alpha=0.7, edgecolors="black", lw=0.4)
_ax_qc.axhline(0, ls="--", lw=0.8, color="grey")
_ax_qc.axvline(0, ls="--", lw=0.8, color="grey")
_ax_qc.axhspan(-_WS_PERM_THR, _WS_PERM_THR, alpha=0.05, color="green")
_ax_qc.axvspan(-_WS_PERM_THR, _WS_PERM_THR, alpha=0.05, color="green")
_ax_qc.set_xlabel("Median permuted log2FE (expected ≈ 0)")
_ax_qc.set_ylabel("Anti-signal log2FE (expected ≈ 0)")
_ax_qc.set_title("Within-sample enrichment sanity checks")
plt.tight_layout()
_fig_qc_ws.savefig(snakemake.output.within_sample_qc_png, dpi=150)
plt.close(_fig_qc_ws)

print(
    f"[Module 4] Within-sample enrichment complete — "
    f"{int(ws_df_out['enriched'].fillna(False).sum())} "
    f"(sample × cell_type) pairs called enriched."
)

print("[Module 4] Done.")
