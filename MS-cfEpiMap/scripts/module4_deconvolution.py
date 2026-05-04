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

# ─── Step 5: Anchor-normalise and compute Poisson z-scores ───────────────────
ctrl_samples = meta.loc[meta["group"] == "Ctrl", "sample_id"].tolist()

rows = []
for ct, sample_cnts in sig_counts.items():
    ctrl_norms = [
        sample_cnts[s] * sf.get(s, 1.0)
        for s in ctrl_samples if s in sample_cnts
    ]
    ctrl_mean = np.mean(ctrl_norms) if ctrl_norms else 1.0

    for sid, cnt in sample_cnts.items():
        norm_cnt = cnt * sf.get(sid, 1.0)
        z = (norm_cnt - ctrl_mean) / np.sqrt(max(ctrl_mean, 1))
        rows.append({
            "sample_id":  sid,
            "cell_type":  ct,
            "raw_count":  cnt,
            "norm_count": norm_cnt,
            "zscore":     z,
        })

scores_df = pd.DataFrame(rows).merge(meta[["sample_id", "group"]], on="sample_id")

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

with pdf_backend.PdfPages(snakemake.output.violin) as pp:
    for ct in cell_types_list:
        ct_df = scores_df[scores_df["cell_type"] == ct].copy()
        if ct_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 5))

        # violin + stripplot
        plot_df = ct_df[ct_df["group"].isin(present_groups)].copy()
        plot_df["group"] = pd.Categorical(plot_df["group"], categories=present_groups)
        plot_df = plot_df.sort_values("group")

        group_data = [plot_df.loc[plot_df["group"] == g, "zscore"].dropna().values
                      for g in present_groups]
        colors = [group_palette.get(g, "grey") for g in present_groups]

        parts = ax.violinplot(group_data, positions=range(len(present_groups)),
                              showmedians=True, showextrema=True)
        for i, (pc, col) in enumerate(zip(parts["bodies"], colors)):
            pc.set_facecolor(col)
            pc.set_alpha(0.6)
        parts["cmedians"].set_color("black")
        parts["cmins"].set_color("black")
        parts["cmaxes"].set_color("black")
        parts["cbars"].set_color("black")

        # jitter
        rng = np.random.default_rng(42)
        for i, (grp, col) in enumerate(zip(present_groups, colors)):
            vals = plot_df.loc[plot_df["group"] == grp, "zscore"].dropna()
            jx   = rng.uniform(-0.08, 0.08, size=len(vals)) + i
            ax.scatter(jx, vals, color=col, s=20, zorder=3, alpha=0.8, edgecolors="white", lw=0.4)

        ax.set_xticks(range(len(present_groups)))
        ax.set_xticklabels([g.replace("MS-Rituximab-", "RTX-") for g in present_groups],
                           rotation=20, ha="right", fontsize=9)
        ax.axhline(0, ls="--", lw=0.8, color="grey")
        ax.set_ylabel("Poisson z-score vs Ctrl")
        ax.set_title(ct)

        # Kruskal-Wallis overall test
        valid = [d for d in group_data if len(d) >= 2]
        if len(valid) >= 2:
            kw_stat, kw_p = stats.kruskal(*valid)
            ax.text(0.98, 0.98, f"KW p={kw_p:.3g}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, color="black")

            # pairwise Dunn annotations on significant pairs only
            posthoc = dunn_posthoc(plot_df, "group", "zscore", present_groups)
            sig_pairs = {(g1, g2): p for (g1, g2), p in posthoc.items() if p < 0.05}
            y_top = plot_df["zscore"].max() if not plot_df.empty else 3
            step  = (y_top - ax.get_ylim()[0]) * 0.12
            for lvl, ((g1, g2), p) in enumerate(sig_pairs.items()):
                x1, x2 = present_groups.index(g1), present_groups.index(g2)
                y = y_top + step * (lvl + 1)
                ax.plot([x1, x1, x2, x2], [y - step*0.2, y, y, y - step*0.2],
                        lw=1, color="black")
                ax.text((x1 + x2) / 2, y, sig_label(p), ha="center", va="bottom", fontsize=9)

        plt.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

    # Summary page: composite indices CDI / NII / IAI
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
            pc.set_facecolor(col)
            pc.set_alpha(0.6)
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
        pp.savefig(fig)
        plt.close(fig)

print("[Module 4] Done.")
