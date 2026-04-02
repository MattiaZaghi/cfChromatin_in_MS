#!/usr/bin/env python3
# path: SNPFingerprintPlot.py
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from decimal import Decimal, getcontext, InvalidOperation

def load_pvalue_matrix(file_path, p_min_str=None, logp_max=None):
    getcontext().prec = 100  # alta precisão

    with open(file_path, 'r') as file:
        lines = file.readlines()

    if lines and lines[0].startswith("#"):
        lines = lines[1:]

    with open(file_path, 'w') as file:
        file.writelines(lines)

    df = pd.read_csv(file_path, sep="\t", index_col=0, dtype=str)

    for suf in [
        "_unique.sorted.dedup.bam",
        ".dac_filtered.dedup.unique.sorted.bam",
    ]:
        df.index = df.index.str.replace(suf, "", regex=False)
        df.columns = df.columns.str.replace(suf, "", regex=False)

    p_floor = None
    if p_min_str is not None:
        try:
            p_floor = Decimal(str(p_min_str).strip())
            if p_floor <= 0 or p_floor >= 1:
                raise ValueError("--p-min must be in (0,1)")
        except (InvalidOperation, ValueError) as e:
            raise ValueError(f"Invalid --p-min: {p_min_str} ({e})")

    def precise_log10_conversion(val_str):
        try:
            s = str(val_str).strip()
            d = Decimal(s)
            if d <= 0:
                return 1000.0
            if p_floor is not None and d < p_floor:
                d = p_floor
            log10_val = d.log10()
            return float(-log10_val)
        except Exception:
            return 1000.0

    df_float = df.applymap(precise_log10_conversion)

    if logp_max is not None:
        df_float = df_float.clip(upper=float(logp_max))

    return df, df_float

def plot_pvalue_heatmap_sorted_dend(df_original, output_path="heatmap.jpg", logp_max=None, vmax_quantile=None, cmap="viridis"):
    if not np.isfinite(df_original.values).all():
        finite = np.asarray(df_original.values)[np.isfinite(df_original.values)]
        fillv = float(np.nanmax(finite)) if finite.size else 0.0
        df_original = df_original.replace([np.inf, -np.inf], np.nan).fillna(fillv)

    row_dist = pdist(df_original, metric='correlation')
    col_dist = pdist(df_original.T, metric='correlation')

    if not np.isfinite(row_dist).all() or not np.isfinite(col_dist).all():
        print("Error: Distance calculation failed, cannot perform clustering")
        return

    row_linkage = linkage(row_dist, method='average')
    col_linkage = linkage(col_dist, method='average')
    size = max(df_original.shape[0]/2,13)
    figsize=(size, size)

    vmin = 0.0
    finite_vals = np.asarray(df_original.values)[np.isfinite(df_original.values)]
    if finite_vals.size == 0:
        vmax = 1.0
    elif logp_max is not None:
        vmax = float(logp_max)
    elif vmax_quantile is not None:
        q = float(vmax_quantile)
        if 0 < q < 1:
            vmax = float(np.quantile(finite_vals, q))
        else:
            vmax = float(np.nanmax(finite_vals))
    else:
        vmax = float(np.nanmax(finite_vals))
    if vmax <= vmin:
        vmax = vmin + 1e-6

    g = sns.clustermap(
        df_original,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        annot=False,
        fmt=".1e",
        cbar_kws={'label': '-log10 (p-value)'},
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        figsize=figsize
    )

    cbar = g.ax_heatmap.collections[0].colorbar
    cbar.set_label('-log10(p-value)', rotation=270, labelpad=15)

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to {output_path}")
    print(f"Color scale limits: vmin={vmin:.3f}, vmax={vmax:.3f}, cmap={cmap}")

def main():
    parser = argparse.ArgumentParser(description="Generate a clustered heatmap from a p-value matrix")
    parser.add_argument("-i", "--input", type=str, default="pval_out.txt", help="Input file path")
    parser.add_argument("-o", "--output", type=str, default="heatmap.jpg", help="Output image file path")
    parser.add_argument("--p-min", type=str, default=None, help="Floor for p-values before transform, e.g. 1e-50")
    parser.add_argument("--logp-max", type=float, default=None, help="Ceiling for -log10(p), e.g. 50")
    parser.add_argument("--vmax-quantile", type=float, default=None, help="Quantile (0-1) for vmax, e.g. 0.99")
    parser.add_argument("--cmap", type=str, default="viridis", help="Colormap name (e.g. viridis, mako, rocket, coolwarm)")
    args = parser.parse_args()

    _, df_float = load_pvalue_matrix(args.input, p_min_str=args.p_min, logp_max=args.logp_max)
    plot_pvalue_heatmap_sorted_dend(df_float, output_path=args.output, logp_max=args.logp_max, vmax_quantile=args.vmax_quantile, cmap=args.cmap)

if __name__ == "__main__":
    main()