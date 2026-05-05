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

# === within-sample enrichment (group-agnostic) ===
#
# Implements per-sample, per-cell-type enrichment that does NOT depend on the
# control cohort.  Background = constitutive anchor regions (≥80 % of reference
# cell types active by construction in Module 1/2 — biologically invariant,
# ideal within-sample null).  For each (sample, cell_type):
#   1. Match n=10 anchor background regions per signature region by length±GC bin.
#   2. Compute log2 fold-enrichment from median densities (reads/kb).
#   3. Fit within-sample NB GLM with length offset → beta_nb, p_nb.
#   4. BH-correct p_nb across cell types within each sample.
#   5. Normalise positive excess densities to fractional contributions per sample.
# Sanity checks:
#   (a) Permuted-signature null (n=200, from anchor regions):
#       |median permuted log2FE| should be < 0.2 per sample.
#   (b) Negative-control anti-signature (random non-RRE windows):
#       log2FE should hover near 0 across the cohort.

# ── WS-1: Load constitutive anchor regions as per-region DataFrame ────────────
# Same region_id format as Module 2 ("chr:start-end") so rows align to anchor_mat.

_anc_df = pd.read_csv(
    snakemake.input.anchors, sep="\t", header=None,
    usecols=[0, 1, 2], names=["chr", "start", "end"],
)
_anc_df["region_id"] = (
    _anc_df["chr"].astype(str) + ":" +
    _anc_df["start"].astype(str) + "-" +
    _anc_df["end"].astype(str)
)
_anc_df["length_bp"] = _anc_df["end"] - _anc_df["start"]

_N_ANCHORS_WS = len(_anc_df)
if _N_ANCHORS_WS < 2000:
    raise RuntimeError(
        f"[Module 4] ABORT: constitutive anchor set has only {_N_ANCHORS_WS} regions "
        "(need ≥2000 for reliable background sampling). "
        "Increase min_anchor_celltypes or add more reference cell types."
    )
print(f"[Module 4] Within-sample: {_N_ANCHORS_WS} constitutive anchor regions.")

# Align to anchor_mat index (only keep regions present in both)
_common_anc_mask = _anc_df["region_id"].isin(anchor_mat.index)
_anc_df = _anc_df[_common_anc_mask].reset_index(drop=True)
print(f"[Module 4] {len(_anc_df)} anchor regions aligned to anchor_mat.")

# ── WS-2: GC content via bedtools nuc (optional) ─────────────────────────────
# If genome FASTA is available, compute pct_gc for every region used in
# background matching.  Otherwise, fall back to length-only matching.
# bedtools nuc 4-col BED output: col0-3 = original, col4 = pct_at, col5 = pct_gc.

_genome_fasta_ws = getattr(snakemake.params, "genome_fasta", "")
_gc_available_ws = False
_gc_map_ws: dict[str, float] = {}


def _bedtools_nuc_gc(bed_df: pd.DataFrame, fasta: str) -> dict[str, float]:
    """Return {region_id: gc_fraction} for rows in bed_df using bedtools nuc."""
    tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    for _, r in bed_df.iterrows():
        tmp.write(f"{r['chr']}\t{int(r['start'])}\t{int(r['end'])}\t{r['region_id']}\n")
    tmp.close()
    res = subprocess.run(
        f"bedtools nuc -fi {fasta} -bed {tmp.name}",
        shell=True, capture_output=True, text=True,
    )
    os.remove(tmp.name)
    gc: dict[str, float] = {}
    for line in res.stdout.splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) > 5:
            try:
                gc[parts[3]] = float(parts[5])
            except (ValueError, IndexError):
                pass
    return gc


if _genome_fasta_ws and Path(_genome_fasta_ws).exists():
    print("[Module 4] Computing GC content via bedtools nuc ...")
    # Collect all unique signature + anchor regions for a single nuc call
    _sig_bed_rows = []
    for ct_ws, sig_df_ws in ct_regions.items():
        if sig_df_ws.empty:
            continue
        for _, row_ws in sig_df_ws.iterrows():
            rid_ws = f"{row_ws[0]}:{int(row_ws[1])}-{int(row_ws[2])}"
            _sig_bed_rows.append({
                "chr": str(row_ws[0]), "start": int(row_ws[1]),
                "end": int(row_ws[2]), "region_id": rid_ws,
            })
    _all_for_gc = pd.concat([
        _anc_df[["chr", "start", "end", "region_id"]],
        pd.DataFrame(_sig_bed_rows).drop_duplicates("region_id"),
    ], ignore_index=True).drop_duplicates("region_id")
    _gc_map_ws = _bedtools_nuc_gc(_all_for_gc, _genome_fasta_ws)
    _gc_available_ws = len(_gc_map_ws) > 0
    print(f"[Module 4] GC computed for {len(_gc_map_ws)} regions.")
else:
    print("[Module 4] No genome FASTA; using length-only background matching.")

_anc_df["gc"] = _anc_df["region_id"].map(_gc_map_ws).fillna(0.5)

# ── WS-3: Per-cell-type region metadata and unique IDs ───────────────────────
# col4 of the per-region counting BED = "{ct}|{chrom}:{start}-{end}"
# This unique ID lets us recover cell type and region after parallel counting.

ws_sig_reg_ids:  dict[str, list]        = {}
ws_sig_lengths:  dict[str, np.ndarray]  = {}
ws_sig_gc_arr:   dict[str, np.ndarray]  = {}
_per_region_rows: list[str] = []

for ct_ws, sig_df_ws in ct_regions.items():
    if sig_df_ws.empty:
        ws_sig_reg_ids[ct_ws]  = []
        ws_sig_lengths[ct_ws]  = np.array([])
        ws_sig_gc_arr[ct_ws]   = np.array([])
        continue
    rids_ct: list[str] = []
    for _, row_ws in sig_df_ws.iterrows():
        rid = f"{ct_ws}|{row_ws[0]}:{int(row_ws[1])}-{int(row_ws[2])}"
        rids_ct.append(rid)
        _per_region_rows.append(
            f"{row_ws[0]}\t{int(row_ws[1])}\t{int(row_ws[2])}\t{rid}"
        )
    ws_sig_reg_ids[ct_ws] = rids_ct
    ws_sig_lengths[ct_ws] = (sig_df_ws[2] - sig_df_ws[1]).values.astype(float) / 1000.0  # kb
    if _gc_available_ws:
        _plain = [f"{row_ws[0]}:{int(row_ws[1])}-{int(row_ws[2])}"
                  for _, row_ws in sig_df_ws.iterrows()]
        ws_sig_gc_arr[ct_ws] = np.array([_gc_map_ws.get(p, 0.5) for p in _plain])
    else:
        ws_sig_gc_arr[ct_ws] = np.full(len(sig_df_ws), 0.5)

# ── WS-4: Anti-signature (negative-control cell type) ────────────────────────
# Random 500-bp genomic windows NOT overlapping the RRE universe.
# Expected log2FE ≈ 0 across samples (no enrichment for arbitrary windows).

_ANTISIG_N   = 500
_ANTISIG_LEN = 500
_chrom_sizes_ws = getattr(snakemake.params, "chrom_sizes", "")
_antisig_ids: list[str] = []
_antisig_len_kb = float(_ANTISIG_LEN) / 1000.0

if _chrom_sizes_ws and Path(_chrom_sizes_ws).exists():
    try:
        _comp_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
        _rand_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
        subprocess.run(
            f"bedtools complement -i {snakemake.input.rre_universe} "
            f"-g {_chrom_sizes_ws} | awk '$3-$2>={_ANTISIG_LEN}' > {_comp_tmp}",
            shell=True, check=True,
        )
        subprocess.run(
            f"bedtools random -g {_chrom_sizes_ws} -n {_ANTISIG_N * 5} "
            f"-l {_ANTISIG_LEN} -seed 123 | "
            f"bedtools intersect -a - -b {_comp_tmp} -u | "
            f"head -{_ANTISIG_N} > {_rand_tmp}",
            shell=True, check=True,
        )
        os.remove(_comp_tmp)
        with open(_rand_tmp) as _fh_as:
            for _line_as in _fh_as:
                _p = _line_as.strip().split("\t")
                if len(_p) >= 3:
                    _rid_as = f"antisig|{_p[0]}:{_p[1]}-{_p[2]}"
                    _antisig_ids.append(_rid_as)
                    _per_region_rows.append(
                        f"{_p[0]}\t{_p[1]}\t{_p[2]}\t{_rid_as}"
                    )
        os.remove(_rand_tmp)
        print(f"[Module 4] Anti-signature: {len(_antisig_ids)} non-RRE windows.")
    except Exception as _e_as:
        print(f"[Module 4] WARNING: anti-signature build failed ({_e_as}); QC check skipped.")

# ── WS-5: Per-region counting (new parallel pass, separate from Step 4) ───────

print(f"\n[Module 4] Per-region counting ({len(samples)} samples) ...")
_raw_pr = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
_raw_pr.write("\n".join(_per_region_rows) + "\n")
_raw_pr.close()
_sorted_pr = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
subprocess.run(
    f"bedtools sort -i {_raw_pr.name} > {_sorted_pr}",
    shell=True, check=True,
)
os.remove(_raw_pr.name)

_pr_args_ws = [
    (sid, f"{frags_dir}/{sid}.bed", _sorted_pr)
    for sid in samples
]
_ws_per_region: dict[str, dict[str, int]] = {}
with ProcessPoolExecutor(max_workers=n_cpus) as _pool_pr:
    for _sid_pr, _prdict in _pool_pr.map(_count_per_region_sample, _pr_args_ws):
        _ws_per_region[_sid_pr] = _prdict
        print(f"  per-region done: {_sid_pr}")
os.remove(_sorted_pr)

# ── WS-6: Extract per-cell-type count arrays per sample ─────────────────────

ws_counts_arr: dict[str, dict[str, np.ndarray]] = {}
ws_antisig_cnts: dict[str, np.ndarray] = {}
for sid in samples:
    ws_counts_arr[sid] = {}
    for ct_ws in CELL_TYPES:
        rids_ws = ws_sig_reg_ids.get(ct_ws, [])
        ws_counts_arr[sid][ct_ws] = np.array(
            [_ws_per_region[sid].get(r, 0) for r in rids_ws], dtype=float
        ) if rids_ws else np.array([], dtype=float)
    if _antisig_ids:
        ws_antisig_cnts[sid] = np.array(
            [_ws_per_region[sid].get(r, 0) for r in _antisig_ids], dtype=float
        )

# ── WS-7: Align anchor_mat and build fast numpy lookup ───────────────────────

_anc_rids_ws   = _anc_df["region_id"].values
_anc_lens_kb   = _anc_df["length_bp"].values.astype(float) / 1000.0
_anc_gc_ws     = _anc_df["gc"].values

_anc_mat_ws = anchor_mat.loc[
    [r for r in _anc_rids_ws if r in anchor_mat.index]
]
_anc_rids_final_ws   = np.array(_anc_mat_ws.index.tolist())
_adf_idx = _anc_df.set_index("region_id")
_anc_lens_final_kb   = _adf_idx.loc[_anc_rids_final_ws, "length_bp"].values.astype(float) / 1000.0
_anc_gc_final_ws     = _adf_idx.loc[_anc_rids_final_ws, "gc"].values

# Convert to numpy array for fast column access: shape (n_anchors, n_samples)
_anc_count_np  = _anc_mat_ws.values.astype(float)
_anc_col_idx   = {s: i for i, s in enumerate(_anc_mat_ws.columns)}
_anc_rid_to_row = {r: i for i, r in enumerate(_anc_rids_final_ws)}

# ── WS-8: Length/GC-matched background sampling ──────────────────────────────
# Bin anchor and signature regions by length decile (and GC decile if available).
# Compute decile edges from the anchor distribution, then assign both.
# For each signature region, sample N_PER_SIG_REGION=10 anchor regions from the
# same (len_bin, gc_bin) cell with replacement.  The result is matched_bg[ct]:
# a flat list of anchor region_ids (length = n_sig_regions × 10).

_WS_N_BINS = 10
_WS_N_PER  = 10
_rng_ws_bg = np.random.default_rng(42)


def _assign_bins(ref_vals: np.ndarray, query_vals: np.ndarray,
                 n_bins: int = 10) -> tuple[np.ndarray, np.ndarray]:
    """Return (ref_bins, query_bins) using percentile edges from ref_vals."""
    qs    = np.linspace(0, 100, n_bins + 1)
    edges = np.percentile(ref_vals, qs)
    edges[-1] += 1e-9  # make right edge inclusive
    r_bins = np.searchsorted(edges[1:], ref_vals)
    q_bins = np.searchsorted(edges[1:], query_vals)
    return r_bins, q_bins


_anc_len_bin_ws, _ = _assign_bins(_anc_lens_final_kb, _anc_lens_final_kb)
if _gc_available_ws:
    _anc_gc_bin_ws, _ = _assign_bins(_anc_gc_final_ws, _anc_gc_final_ws)
    _anc_bin_key_ws = _anc_len_bin_ws * _WS_N_BINS + _anc_gc_bin_ws
else:
    _anc_bin_key_ws = _anc_len_bin_ws

_bin_to_row_idx: dict[int, list[int]] = defaultdict(list)
for _i_anc, _bk_anc in enumerate(_anc_bin_key_ws):
    _bin_to_row_idx[int(_bk_anc)].append(_i_anc)
_all_anc_rows = list(range(len(_anc_rids_final_ws)))

matched_bg_ws: dict[str, list[int]] = {}  # ct → list of anchor row indices
for ct_ws in CELL_TYPES:
    _lens_ct = ws_sig_lengths.get(ct_ws, np.array([]))  # already in kb
    _gc_ct   = ws_sig_gc_arr.get(ct_ws, np.array([]))
    if len(_lens_ct) == 0:
        matched_bg_ws[ct_ws] = []
        continue
    _, _sig_len_b = _assign_bins(_anc_lens_final_kb, _lens_ct)
    if _gc_available_ws:
        _, _sig_gc_b = _assign_bins(_anc_gc_final_ws, _gc_ct)
        _sig_bk = _sig_len_b * _WS_N_BINS + _sig_gc_b
    else:
        _sig_bk = _sig_len_b
    _bg_row_indices: list[int] = []
    for _bk_s in _sig_bk:
        _pool = _bin_to_row_idx.get(int(_bk_s), _all_anc_rows)
        _chosen = _rng_ws_bg.choice(_pool, size=_WS_N_PER, replace=True)
        _bg_row_indices.extend(_chosen.tolist())
    matched_bg_ws[ct_ws] = _bg_row_indices

# ── WS-9: Per-sample NB GLM and log2FE ───────────────────────────────────────

_WS_EPS     = 1.0       # pseudocount for log2FE to avoid log(0)
_WS_MIN_SIG = 50        # skip NB GLM below this region count


def _fit_nb_glm_ws(
    sig_cnt: np.ndarray, bg_cnt: np.ndarray,
    sig_len_kb: np.ndarray, bg_len_kb: np.ndarray,
) -> tuple[float, float]:
    """
    NB GLM with log(length_kb) offset.  is_sig=1 for signature, 0 for background.
    Returns (beta_nb, p_nb) for the is_sig coefficient, or (nan, nan) on failure.
    """
    counts  = np.r_[sig_cnt, bg_cnt]
    lengths = np.maximum(np.r_[sig_len_kb, bg_len_kb], 1e-3)
    is_sig  = np.r_[np.ones(len(sig_cnt)), np.zeros(len(bg_cnt))]
    X = sm.add_constant(is_sig, has_constant="add")
    try:
        fit = sm.GLM(
            counts, X,
            family=sm.families.NegativeBinomial(),
            offset=np.log(lengths),
        ).fit(disp=False, maxiter=100)
        return float(fit.params[1]), float(fit.pvalues[1])
    except Exception:
        return np.nan, np.nan


ws_rows_out: list[dict] = []
for sid in samples:
    _col_ws = _anc_col_idx.get(sid, -1)
    for ct_ws in CELL_TYPES:
        sig_cnt = ws_counts_arr[sid].get(ct_ws, np.array([]))
        sig_len = ws_sig_lengths.get(ct_ws, np.array([]))  # kb
        bg_rows = matched_bg_ws.get(ct_ws, [])

        if len(sig_cnt) == 0 or len(bg_rows) == 0:
            ws_rows_out.append({
                "sample": sid, "cell_type": ct_ws,
                "log2FE": np.nan, "beta_nb": np.nan, "p_nb": np.nan,
                "excess": np.nan, "low_n": True,
            })
            continue

        _bg_rows_arr = np.array(bg_rows, dtype=int)
        bg_cnt = (
            _anc_count_np[_bg_rows_arr, _col_ws]
            if _col_ws >= 0 else np.zeros(len(_bg_rows_arr))
        )
        bg_len = _anc_lens_final_kb[_bg_rows_arr]

        sig_dens = sig_cnt / np.maximum(sig_len, 1e-3)
        bg_dens  = bg_cnt  / np.maximum(bg_len,  1e-3)
        med_sig  = float(np.median(sig_dens))
        med_bg   = float(np.median(bg_dens))
        log2fe   = float(np.log2((med_sig + _WS_EPS) / (med_bg + _WS_EPS)))
        excess   = max(0.0, med_sig - med_bg)

        low_n = len(sig_cnt) < _WS_MIN_SIG
        if not low_n:
            beta_nb, p_nb = _fit_nb_glm_ws(sig_cnt, bg_cnt, sig_len, bg_len)
        else:
            beta_nb, p_nb = np.nan, np.nan

        ws_rows_out.append({
            "sample": sid, "cell_type": ct_ws,
            "log2FE": log2fe, "beta_nb": beta_nb, "p_nb": p_nb,
            "excess": excess, "low_n": low_n,
        })

ws_df_out = pd.DataFrame(ws_rows_out)


def _ws_bh_fdr(grp: pd.DataFrame) -> pd.DataFrame:
    pvals = grp["p_nb"].values
    valid = ~np.isnan(pvals)
    q = np.full(len(pvals), np.nan)
    if valid.sum() > 1:
        q[valid] = _bh_correct(pvals[valid])
    elif valid.sum() == 1:
        q[valid] = pvals[valid]
    grp = grp.copy()
    grp["q_nb"] = q
    grp["enriched"] = grp["q_nb"] < 0.05
    return grp


ws_df_out = ws_df_out.groupby("sample", group_keys=False).apply(_ws_bh_fdr)


def _ws_contrib(grp: pd.DataFrame) -> pd.DataFrame:
    exc   = grp["excess"].clip(lower=0).fillna(0)
    total = exc.sum()
    grp = grp.copy()
    grp["contribution"] = (exc / total) if total > 0 else 0.0
    return grp


ws_df_out = ws_df_out.groupby("sample", group_keys=False).apply(_ws_contrib)

# ── WS-10: Sanity checks ──────────────────────────────────────────────────────

_N_PERM_WS = 200
_rng_perm_ws = np.random.default_rng(99)
_n_anc_ws = len(_anc_rids_final_ws)
qc_ws_rows: list[dict] = []

for sid in samples:
    _col_ws = _anc_col_idx.get(sid, -1)
    _ac_ws = (
        _anc_count_np[:, _col_ws]
        if _col_ws >= 0 else np.zeros(_n_anc_ws)
    )
    _bg_med_ws = float(np.median(_ac_ws / np.maximum(_anc_lens_final_kb, 1e-3)))

    # Permuted null: sample size-matched sets from anchor regions
    _perm_log2fes: list[float] = []
    _perm_size = max(50, _n_anc_ws // 10)
    for _ in range(_N_PERM_WS):
        _idx_p = _rng_perm_ws.integers(0, _n_anc_ws, size=_perm_size)
        _p_med = float(np.median(_ac_ws[_idx_p] / np.maximum(_anc_lens_final_kb[_idx_p], 1e-3)))
        _perm_log2fes.append(np.log2((_p_med + _WS_EPS) / (_bg_med_ws + _WS_EPS)))
    _med_perm = float(np.median(_perm_log2fes))

    # Anti-signature
    if _antisig_ids and sid in ws_antisig_cnts:
        _as_vec = ws_antisig_cnts[sid]
        _as_len = np.full(len(_as_vec), _antisig_len_kb)
        _as_med = float(np.median(_as_vec / np.maximum(_as_len, 1e-3)))
        _neg_ctrl_log2fe = float(np.log2((_as_med + _WS_EPS) / (_bg_med_ws + _WS_EPS)))
    else:
        _neg_ctrl_log2fe = np.nan

    _n_enr = int(ws_df_out.loc[ws_df_out["sample"] == sid, "enriched"].fillna(False).sum())
    qc_ws_rows.append({
        "sample": sid,
        "median_perm_log2FE":    _med_perm,
        "neg_ctrl_log2FE":       _neg_ctrl_log2fe,
        "n_signatures_enriched": _n_enr,
    })

qc_ws_df = pd.DataFrame(qc_ws_rows)

_suspect_ws = qc_ws_df[qc_ws_df["median_perm_log2FE"].abs() > 0.2]
if len(_suspect_ws) > 0:
    print(
        f"[Module 4] WARNING: {len(_suspect_ws)} sample(s) have "
        f"|median_perm_log2FE| > 0.2 — possible systematic bias:\n"
        + _suspect_ws[["sample", "median_perm_log2FE"]].to_string(index=False)
    )
if len(_suspect_ws) == len(samples) and len(samples) > 0:
    raise RuntimeError(
        "[Module 4] ABORT: permuted-signature null fails for every sample. "
        "Within-sample enrichment is unreliable. "
        "Check anchor regions and signature BEDs."
    )

# ── WS-11: Write TSV outputs ──────────────────────────────────────────────────

_ws_out_dir = Path(snakemake.output.within_sample_scores).parent
_ws_out_dir.mkdir(parents=True, exist_ok=True)

ws_df_out[[
    "sample", "cell_type", "log2FE", "beta_nb", "p_nb", "q_nb",
    "excess", "contribution", "enriched", "low_n",
]].to_csv(snakemake.output.within_sample_scores, sep="\t", index=False)

qc_ws_df.to_csv(snakemake.output.within_sample_qc, sep="\t", index=False)

# ── WS-12: Figures ────────────────────────────────────────────────────────────

_ct_ws_order = list(CELL_TYPES.keys())

# Fig log2FE: group-agnostic boxplot (page 1) + group-coloured scatter (page 2)
_box_data = [
    ws_df_out.loc[ws_df_out["cell_type"] == ct, "log2FE"].dropna().values
    for ct in _ct_ws_order
]
fig_log2fe_ws, ax_log2fe = plt.subplots(figsize=(max(8, len(_ct_ws_order) * 0.6), 5))
ax_log2fe.boxplot(
    _box_data, labels=_ct_ws_order,
    patch_artist=True, notch=False, sym=".",
    medianprops={"color": "black", "linewidth": 1.5},
    boxprops={"facecolor": "#AECDE8", "alpha": 0.8},
)
ax_log2fe.axhline(0, ls="--", lw=0.8, color="grey")
ax_log2fe.set_ylabel("log2 fold enrichment over anchor baseline")
ax_log2fe.set_title("Within-sample tissue enrichment (group-agnostic)")
ax_log2fe.set_xticklabels(_ct_ws_order, rotation=45, ha="right", fontsize=8)
plt.tight_layout()

# Group-coloured variant (second page — downstream reference)
_ws_grp_df = ws_df_out.merge(
    meta[["sample_id", "group"]], left_on="sample", right_on="sample_id", how="left"
)
_n_ct = len(_ct_ws_order)
fig_grp_ws, axes_grp_ws = plt.subplots(1, _n_ct, figsize=(max(12, _n_ct * 1.1), 5), sharey=True)
if _n_ct == 1:
    axes_grp_ws = [axes_grp_ws]
_rng_jitter_ws = np.random.default_rng(7)
for _ax_g, _ct_g in zip(axes_grp_ws, _ct_ws_order):
    _ct_sub_g = _ws_grp_df[_ws_grp_df["cell_type"] == _ct_g]
    for _grp_g, _col_g in group_palette.items():
        _vals_g = _ct_sub_g.loc[_ct_sub_g["group"] == _grp_g, "log2FE"].dropna()
        if len(_vals_g) > 0:
            _jx_g = _rng_jitter_ws.uniform(-0.15, 0.15, size=len(_vals_g))
            _ax_g.scatter(_jx_g, _vals_g, color=_col_g, s=12, alpha=0.7, label=_grp_g)
    _ax_g.axhline(0, ls="--", lw=0.6, color="grey")
    _ax_g.set_title(_ct_g, fontsize=7, rotation=45, ha="right")
    _ax_g.set_xticks([])
_handles_ws = [
    plt.Line2D([0], [0], marker="o", color=c, lw=0, label=g)
    for g, c in group_palette.items()
]
fig_grp_ws.legend(handles=_handles_ws, loc="upper right", fontsize=7)
fig_grp_ws.suptitle("Within-sample log2FE by group (downstream reference only)")
plt.tight_layout()

with pdf_backend.PdfPages(snakemake.output.within_sample_log2fe) as _pp_ws:
    _pp_ws.savefig(fig_log2fe_ws)  # page 1: group-agnostic
    _pp_ws.savefig(fig_grp_ws)     # page 2: group-coloured variant
plt.close(fig_log2fe_ws)
plt.close(fig_grp_ws)

# Fig contribution: per-sample stacked bars, alphabetical sample order
_contrib_piv = ws_df_out.pivot(
    index="sample", columns="cell_type", values="contribution"
).fillna(0).sort_index()
_ws_ct_colors = sns.color_palette("tab20", n_colors=len(_ct_ws_order))
_ws_ct_cmap   = dict(zip(_ct_ws_order, _ws_ct_colors))

fig_contrib_ws, ax_contrib = plt.subplots(
    figsize=(max(10, len(_contrib_piv) * 0.15), 5)
)
_bot = np.zeros(len(_contrib_piv))
for _ct_cb in _ct_ws_order:
    if _ct_cb not in _contrib_piv.columns:
        continue
    _vals_cb = _contrib_piv[_ct_cb].values
    ax_contrib.bar(
        range(len(_contrib_piv)), _vals_cb, bottom=_bot,
        color=_ws_ct_cmap[_ct_cb], label=_ct_cb, width=1.0, edgecolor="none",
    )
    _bot += _vals_cb
ax_contrib.set_xlim(-0.5, len(_contrib_piv) - 0.5)
ax_contrib.set_ylim(0, 1)
ax_contrib.set_ylabel("Fractional contribution")
ax_contrib.set_title("Per-sample cell-type contributions (alphabetical sample order)")
ax_contrib.set_xticks([])
ax_contrib.set_xlabel("Samples (alphabetical order, no group sorting)")
ax_contrib.legend(loc="upper right", fontsize=6, ncol=2)
plt.tight_layout()
with pdf_backend.PdfPages(snakemake.output.within_sample_contribution) as _pp_cb:
    _pp_cb.savefig(fig_contrib_ws)
plt.close(fig_contrib_ws)

# QC PNG: scatter of permuted null vs anti-signature log2FE per sample
fig_qc_ws, ax_qc = plt.subplots(figsize=(6, 5))
ax_qc.scatter(
    qc_ws_df["median_perm_log2FE"], qc_ws_df["neg_ctrl_log2FE"],
    s=30, alpha=0.7, edgecolors="black", lw=0.4,
)
ax_qc.axhline(0, ls="--", lw=0.8, color="grey")
ax_qc.axvline(0, ls="--", lw=0.8, color="grey")
ax_qc.axhspan(-0.2, 0.2, alpha=0.05, color="green")
ax_qc.axvspan(-0.2, 0.2, alpha=0.05, color="green")
ax_qc.set_xlabel("Median permuted log2FE (expected ≈ 0)")
ax_qc.set_ylabel("Anti-signature log2FE (expected ≈ 0)")
ax_qc.set_title("Within-sample enrichment sanity checks")
plt.tight_layout()
fig_qc_ws.savefig(snakemake.output.within_sample_qc_png, dpi=150)
plt.close(fig_qc_ws)

print(
    f"[Module 4] Within-sample enrichment complete — "
    f"{int(ws_df_out['enriched'].fillna(False).sum())} "
    f"(sample × cell_type) pairs called enriched."
)

print("[Module 4] Done.")
