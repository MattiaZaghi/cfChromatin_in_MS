"""
Module 8 — Transcription Factor Activity (AUC-based, optimized).

KEY OPTIMIZATION vs previous version:
  Old: N_TFs × N_samples × 2 serial bedtools calls, each re-sorting the
       fragment BED (e.g. 14 TFs × 56 samples × 2 = 1568 sorts).
  New: Combine all TF sites + anchor regions into one labeled BED.
       Sort midpoints ONCE per sample, one intersect covers everything.
       Parallelized across samples → 56 parallel jobs total.
       Typical runtime: 5-15 min instead of hours.
"""

import glob
import os
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend


_ANCHOR_LABEL = "__ANCHOR__"


# ─── Parallel worker ──────────────────────────────────────────────────────────

def _count_one_sample(args):
    """
    Sort midpoints ONCE, intersect with combined TF+anchor BED.
    Returns (sid, {label: np.array of per-region counts}).
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

    label_counts: dict[str, list] = defaultdict(list)
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 5:
            label_counts[parts[3]].append(int(parts[-1]))

    return sid, {k: np.array(v, dtype=float) for k, v in label_counts.items()}


# ─── Step 1: Load inputs ──────────────────────────────────────────────────────
samples      = snakemake.params.samples
frags_dir    = snakemake.params.frags_dir
min_sites    = int(snakemake.params.min_sites)
tf_sites_dir = str(snakemake.input.tf_sites_dir)
anchors_bed  = snakemake.input.anchors
n_cpus       = snakemake.threads

meta  = pd.read_csv(snakemake.input.meta, sep="\t", comment="#")
meta  = meta[meta["sample_id"].isin(samples)]
sf_df = pd.read_csv(snakemake.input.sf, sep="\t")
sf    = sf_df.set_index("sample_id")["constitutive_sf"].to_dict()

# ─── Step 2: Discover TF BED files ───────────────────────────────────────────
tf_beds = {
    Path(p).stem: p
    for p in glob.glob(os.path.join(tf_sites_dir, "*.bed"))
    if sum(1 for _ in open(p)) >= min_sites
}

if not tf_beds:
    print(
        f"[Module 8] WARNING: No TF BED files with >={min_sites} sites found "
        f"in {tf_sites_dir}. Writing placeholder outputs."
    )
    Path(snakemake.output.auc_matrix).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=["TF"] + list(samples)).to_csv(
        snakemake.output.auc_matrix, sep="\t", index=False
    )
    for out_pdf in [snakemake.output.heatmap, snakemake.output.bcell_plot]:
        with pdf_backend.PdfPages(out_pdf) as pp:
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.text(0.5, 0.5,
                    "No TF BED files found.\nPlace hg19 BEDs in reference/tf_sites/",
                    ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            pp.savefig(fig)
            plt.close(fig)
    raise SystemExit(0)

print(f"[Module 8] {len(tf_beds)} TFs with >={min_sites} sites")

# ─── Step 3: Build combined BED (TF sites + anchors, label in col 4) ─────────
combined_rows = []

for tf_name, tf_bed in tf_beds.items():
    with open(tf_bed) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                combined_rows.append(
                    f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{tf_name}"
                )

with open(anchors_bed) as fh:
    for line in fh:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            combined_rows.append(
                f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{_ANCHOR_LABEL}"
            )

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
print(f"[Module 8] Counting for {len(samples)} samples "
      f"using {n_cpus} parallel workers ...")

args_list = [
    (sid, f"{frags_dir}/{sid}.bed", sorted_combined)
    for sid in samples
]

all_counts: dict[str, dict[str, np.ndarray]] = {}
with ProcessPoolExecutor(max_workers=n_cpus) as pool:
    for sid, label_counts in pool.map(_count_one_sample, args_list):
        all_counts[sid] = label_counts
        print(f"  done: {sid}")

os.remove(sorted_combined)

# ─── Step 5: Compute AUC per TF per sample ────────────────────────────────────
rng = np.random.default_rng(42)
auc_matrix = pd.DataFrame(index=list(tf_beds.keys()), columns=samples, dtype=float)

for sid in samples:
    sample_counts = all_counts.get(sid, {})
    bg_full = sample_counts.get(_ANCHOR_LABEL, np.array([]))

    for tf_name in tf_beds:
        tf_counts = sample_counts.get(tf_name, np.array([]))
        if len(tf_counts) == 0 or len(bg_full) == 0:
            auc_matrix.loc[tf_name, sid] = np.nan
            continue

        # Down-sample background to match TF site count (balanced AUC)
        if len(bg_full) > len(tf_counts):
            bg = rng.choice(bg_full, size=len(tf_counts), replace=False)
        else:
            bg = bg_full

        # Apply scaling factor
        tf_sc = tf_counts * sf.get(sid, 1.0)
        bg_sc = bg        * sf.get(sid, 1.0)

        all_vals   = np.concatenate([tf_sc, bg_sc])
        all_labels = np.concatenate([np.ones(len(tf_sc)), np.zeros(len(bg_sc))])
        try:
            auc_matrix.loc[tf_name, sid] = roc_auc_score(all_labels, all_vals)
        except ValueError:
            auc_matrix.loc[tf_name, sid] = np.nan

# ─── Step 6: Save AUC matrix ─────────────────────────────────────────────────
Path(snakemake.output.auc_matrix).parent.mkdir(parents=True, exist_ok=True)
auc_matrix.to_csv(snakemake.output.auc_matrix, sep="\t")

# ─── Step 7: Heatmap ─────────────────────────────────────────────────────────
groups      = ["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
palette     = {
    "Ctrl": "#4878CF", "NEW": "#6ACC65",
    "MS-Rituximab-Stable": "#D65F5F",
    "MS-Rituximab-Progressive": "#B47CC7",
}
_mi = meta.set_index("sample_id")
group_order  = [g for g in groups if g in meta["group"].unique()]
sample_order = (
    _mi.loc[_mi["group"].isin(group_order), "group"]
    .sort_values().index.tolist()
)
sample_order = [s for s in sample_order if s in auc_matrix.columns]

with pdf_backend.PdfPages(snakemake.output.heatmap) as pp:
    import seaborn as sns
    fig, ax = plt.subplots(
        figsize=(max(8, len(sample_order) * 0.15), max(5, len(tf_beds) * 0.35))
    )
    sns.heatmap(
        auc_matrix[sample_order].astype(float),
        cmap="RdBu_r", center=0.5, vmin=0.3, vmax=0.7,
        ax=ax, yticklabels=True, xticklabels=False,
        cbar_kws={"label": "AUC"},
    )
    ax.set_title("TF Activity (AUC vs constitutive background)")
    ax.set_xlabel("Samples")
    plt.tight_layout()
    pp.savefig(fig)
    plt.close(fig)

# ─── Step 8: B cell TF validation plot ───────────────────────────────────────
bcell_tfs = [tf for tf in ["PAX5", "IRF4", "EBF1"] if tf in auc_matrix.index]

with pdf_backend.PdfPages(snakemake.output.bcell_plot) as pp:
    for tf in bcell_tfs:
        fig, ax = plt.subplots(figsize=(5, 4))
        auc_row = auc_matrix.loc[tf].astype(float)
        merged  = meta[["sample_id", "group"]].copy()
        merged[tf] = merged["sample_id"].map(auc_row)
        for grp, sub in merged.dropna(subset=[tf]).groupby("group"):
            ax.scatter(range(len(sub)), sub[tf],
                       label=grp, color=palette.get(grp, "grey"), alpha=0.8)
        ax.axhline(0.5, ls="--", color="grey")
        ax.set_title(f"{tf} AUC by group")
        ax.set_ylabel("AUC")
        ax.legend(fontsize=7)
        plt.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

print("[Module 8] Done.")
