#!/usr/bin/env python3
"""
plot_background_rates.py
=========================
For every sample, plots the Sadeh local-background-rate model used by
MS-cfEpiMap/scripts/module4_deconvolution.py (WS-F: genome -> chromosome ->
2 Mb regional hierarchy, 95th-percentile-trimmed-mean fitNoise) against
genomic position.

Reuses the exact background-tile universe module4 counts against
(results/deconvolution/background_signal_windows.bed.gz, "bg|" rows) and the
exact hierarchy constants/formula (_WS_REGION, _WS_BG_TRIM, _WS_MIN_BG,
fitNoise), but recomputes per-sample per-tile fragment counts standalone --
module4 does not persist those counts or the fitted rate models anywhere in
results/, only an aggregate fit-fraction diagnostic. Counting is the
expensive step, so results are cached per sample.

One PNG, one row per sample (long & short panels so many stack on one page):
  x = genomic coordinate, chr1..chr22 concatenated, each chromosome's width
      proportional to its hg19 size; red vertical lines at chromosome
      boundaries.
  y = reads/kb.
  grey  = genome-wide rate (single value).
  black = per-chromosome rate (step per chromosome).
  blue  = per-2Mb-region local rate (step per bin).
Runs of >= --gap-min consecutive 2 Mb bins lacking a local estimate (fewer
than --min-bg background tiles in that bin) are left as a visible gap rather
than bridged; shorter runs are linearly interpolated so the line doesn't
needlessly break on a single missing bin. Applies to both the chromosome and
regional lines (whichever actually has missing bins -- the genome line is a
single value and never gaps).

Usage
-----
    python plot_background_rates.py
        [--samples s1,s2,...]      # default: all qc_include==TRUE samples
        [--gap-min 3] [--cpus 4]
        [--out background_rate_plots.png]

Counting ~464k background tiles per sample is comparable in cost to what
module4 already pays per sample (it counts against a combined BED of very
similar size) -- expect it to take a while for the full cohort the first
time; results are cached under --cache-dir so re-running for plot cosmetics
is instant. Each --cpus worker runs its own bedtools/awk pipeline, so
memory scales with --cpus roughly like module4's own per-sample workers.
"""
import argparse
import gzip
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CFCHROM_DIR = Path("/cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS")
DEFAULT_BG_BED = Path(
    "/cfs/klemming/projects/supr/uppstore2017150/anita/cfchrom_frag/"
    "results/deconvolution/background_signal_windows.bed.gz"
)
DEFAULT_META = CFCHROM_DIR / "MS-cfEpiMap/config/sample_metadata.tsv"
DEFAULT_CHROM_SIZES = Path(
    "/cfs/klemming/projects/supr/uppstore2017150/anita/hg19/hg19.chrom.sizes"
)
DEFAULT_OUT = CFCHROM_DIR / "background_rate_plots.png"
DEFAULT_CACHE_DIR = CFCHROM_DIR / "_cache_background_rate_counts"

# Must match module4_deconvolution.py's WS-A constants exactly.
WS_REGION = 2_000_000
WS_BG_TRIM = 0.95
WS_MIN_BG = 20

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def load_bg_tiles(bg_bed_path: Path) -> pd.DataFrame:
    """Load the 'bg|' rows (autosomal off-target tiles) from module4's
    combined counting BED."""
    rows = []
    opener = gzip.open if str(bg_bed_path).endswith(".gz") else open
    with opener(bg_bed_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or not parts[3].startswith("bg|"):
                continue
            chrom, start, end, rid = parts[0], int(parts[1]), int(parts[2]), parts[3]
            if chrom not in AUTOSOMES:
                continue
            rows.append((chrom, start, end, rid))
    df = pd.DataFrame(rows, columns=["chr", "start", "end", "region_id"])
    df["length_kb"] = (df["end"] - df["start"]) / 1000.0
    df["region_2mb"] = df["chr"] + ":" + (df["start"] // WS_REGION).astype(str)
    return df


def count_one_sample(job):
    """Sort midpoints once, intersect once against the sorted bg-tile BED.
    Mirrors module4_deconvolution.py's _count_per_region_sample exactly."""
    sid, frags_bed, sorted_bg_bed = job
    mid_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    try:
        awk = (
            f"zcat -f {frags_bed} | "
            f"awk 'BEGIN{{OFS=\"\\t\"}}"
            f"{{mid=int(($2+$3)/2); print $1,mid,mid+1}}'"
        )
        p1 = subprocess.Popen(awk, shell=True, stdout=subprocess.PIPE)
        p2 = subprocess.run(
            "bedtools sort -i -", shell=True, stdin=p1.stdout, capture_output=True
        )
        p1.stdout.close()
        with open(mid_tmp, "wb") as f:
            f.write(p2.stdout)
        res = subprocess.run(
            f"bedtools intersect -a {sorted_bg_bed} -b {mid_tmp} -c",
            shell=True, capture_output=True, text=True,
        )
    finally:
        if os.path.exists(mid_tmp):
            os.remove(mid_tmp)

    counts = {}
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 5:
            counts[parts[3]] = int(parts[-1])
    return sid, counts


def fit_noise(densities: np.ndarray) -> float:
    """95th-pct-trimmed mean density -- identical to module4's _ws_fit_noise."""
    x = densities[np.isfinite(densities) & (densities >= 0)]
    if len(x) < WS_MIN_BG:
        return float("nan")
    x = x[x <= np.quantile(x, WS_BG_TRIM)]
    return float(np.mean(x)) if len(x) > 0 else float("nan")


def bridge_short_gaps(values: np.ndarray, gap_min: int) -> np.ndarray:
    """Linearly interpolate runs of < gap_min consecutive NaNs (so the line
    still connects the two adjacent real estimates); runs of >= gap_min
    NaNs are left as NaN, which matplotlib renders as a visible gap."""
    values = values.astype(float).copy()
    n = len(values)
    isnan = np.isnan(values)
    i = 0
    while i < n:
        if not isnan[i]:
            i += 1
            continue
        j = i
        while j < n and isnan[j]:
            j += 1
        run_len = j - i
        if run_len < gap_min and i > 0 and j < n:
            left, right = values[i - 1], values[j]
            values[i:j] = np.linspace(left, right, run_len + 2)[1:-1]
        i = j
    return values


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--bg-bed", type=Path, default=DEFAULT_BG_BED)
    ap.add_argument("--sample-metadata", type=Path, default=DEFAULT_META)
    ap.add_argument("--chrom-sizes", type=Path, default=DEFAULT_CHROM_SIZES)
    ap.add_argument(
        "--samples", type=str, default=None,
        help="Comma-separated sample_id subset. Default: all "
             "qc_include==TRUE samples in sample_metadata.tsv.",
    )
    ap.add_argument(
        "--gap-min", type=int, default=3,
        help="Runs of >= this many consecutive 2 Mb bins without a local "
             "estimate are left as a visible gap; shorter runs are bridged.",
    )
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument(
        "--cache-dir", type=Path, default=DEFAULT_CACHE_DIR,
        help="Per-sample bg-tile counts are cached here. Delete this dir "
             "to force a recount.",
    )
    ap.add_argument("--panel-height", type=float, default=1.1, help="Inches per sample row.")
    ap.add_argument("--fig-width", type=float, default=22.0)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    print("[1/5] Loading background tiles ...")
    bg_df = load_bg_tiles(args.bg_bed)
    print(f"  {len(bg_df):,} autosomal background tiles "
          f"({bg_df['chr'].nunique()} chromosomes).")

    bg_bed_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    for r in bg_df.itertuples(index=False):
        bg_bed_tmp.write(f"{r.chr}\t{r.start}\t{r.end}\t{r.region_id}\n")
    bg_bed_tmp.close()
    sorted_bg_bed = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    subprocess.run(f"bedtools sort -i {bg_bed_tmp.name} > {sorted_bg_bed}",
                    shell=True, check=True)
    os.remove(bg_bed_tmp.name)

    print("[2/5] Resolving samples ...")
    meta = pd.read_csv(args.sample_metadata, sep="\t", comment="#")
    if "qc_include" in meta.columns:
        meta = meta[meta["qc_include"].astype(str).str.upper() == "TRUE"]
    if args.samples:
        wanted = {s.strip() for s in args.samples.split(",")}
        meta = meta[meta["sample_id"].isin(wanted)]
    sample_bed_map = dict(zip(meta["sample_id"], meta["fragment_bed"]))
    samples = list(sample_bed_map)
    print(f"  {len(samples)} samples.")

    print(f"[3/5] Counting background-tile midpoints per sample "
          f"({args.cpus} workers; cached in {args.cache_dir}) ...")
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    sample_counts = {}
    to_run = []
    for sid in samples:
        cache_f = args.cache_dir / f"{sid}.tsv.gz"
        if cache_f.exists():
            sample_counts[sid] = pd.read_csv(
                cache_f, sep="\t", index_col=0
            )["count"].to_dict()
        else:
            to_run.append(sid)
    if to_run:
        jobs = [(sid, sample_bed_map[sid], sorted_bg_bed) for sid in to_run]
        with ProcessPoolExecutor(max_workers=args.cpus) as pool:
            for sid, counts in pool.map(count_one_sample, jobs):
                sample_counts[sid] = counts
                pd.Series(counts, name="count").to_csv(
                    args.cache_dir / f"{sid}.tsv.gz", sep="\t"
                )
                print(f"    done: {sid}")
    os.remove(sorted_bg_bed)

    print("[4/5] Fitting per-sample background models + genome axis ...")
    chrom_sizes = {}
    with open(args.chrom_sizes) as fh:
        for line in fh:
            c, sz = line.split()[:2]
            if c in AUTOSOMES:
                chrom_sizes[c] = int(sz)
    chrom_order = [c for c in AUTOSOMES if c in chrom_sizes]
    offset, cum = {}, 0
    for c in chrom_order:
        offset[c] = cum
        cum += chrom_sizes[c]
    genome_len = cum

    bin_chr, bin_start = [], []
    for c in chrom_order:
        n_bins = int(np.ceil(chrom_sizes[c] / WS_REGION))
        for b in range(n_bins):
            bin_chr.append(c)
            bin_start.append(b * WS_REGION)
    bin_chr = np.array(bin_chr)
    bin_start = np.array(bin_start)
    bin_x = np.array([offset[c] + s + WS_REGION / 2 for c, s in zip(bin_chr, bin_start)])
    order_idx = np.argsort(bin_x)
    bin_chr, bin_start, bin_x = bin_chr[order_idx], bin_start[order_idx], bin_x[order_idx]
    region_key = np.array([f"{c}:{s // WS_REGION}" for c, s in zip(bin_chr, bin_start)])

    bg_region_2mb = bg_df["region_2mb"].to_numpy()
    bg_chr = bg_df["chr"].to_numpy()
    bg_len_kb = bg_df["length_kb"].to_numpy()
    bg_rid = bg_df["region_id"].to_numpy()

    models = {}
    for sid in samples:
        cnts = sample_counts.get(sid, {})
        counts_arr = np.array([cnts.get(r, 0) for r in bg_rid], dtype=float)
        density = counts_arr / bg_len_kb

        genome_rate = fit_noise(density)

        chr_rate_map = {c: fit_noise(density[bg_chr == c]) for c in chrom_order}
        chr_rate_per_bin = bridge_short_gaps(
            np.array([chr_rate_map.get(c, np.nan) for c in bin_chr]), args.gap_min
        )

        reg_rate_map = {
            rk: fit_noise(density[bg_region_2mb == rk])
            for rk in np.unique(bg_region_2mb)
        }
        reg_rate_per_bin = bridge_short_gaps(
            np.array([reg_rate_map.get(rk, np.nan) for rk in region_key]), args.gap_min
        )

        models[sid] = {
            "genome": genome_rate,
            "chr_per_bin": chr_rate_per_bin,
            "reg_per_bin": reg_rate_per_bin,
        }

    print("[5/5] Plotting ...")
    n = len(samples)
    fig, axes = plt.subplots(
        n, 1, figsize=(args.fig_width, args.panel_height * n), sharex=True, squeeze=False
    )
    axes = axes[:, 0]

    chrom_bounds = [offset[c] for c in chrom_order] + [genome_len]
    chrom_mid = [
        (chrom_bounds[i] + chrom_bounds[i + 1]) / 2 for i in range(len(chrom_order))
    ]

    for ax, sid in zip(axes, samples):
        m = models[sid]
        ax.axhline(
            m["genome"], color="grey", lw=1.0, zorder=1,
            label="genome" if sid == samples[0] else None,
        )
        ax.plot(
            bin_x, m["chr_per_bin"], color="black", lw=1.0, drawstyle="steps-mid",
            zorder=2, label="chromosome" if sid == samples[0] else None,
        )
        ax.plot(
            bin_x, m["reg_per_bin"], color="#3b6fb0", lw=0.9, drawstyle="steps-mid",
            zorder=3, label="local (2 Mb)" if sid == samples[0] else None,
        )
        for b in chrom_bounds[1:-1]:
            ax.axvline(b, color="red", lw=0.6, zorder=0)
        ax.set_xlim(0, genome_len)
        ax.set_ylabel("reads/kb", fontsize=7, color="black")
        ax.tick_params(axis="y", labelsize=6, colors="black")
        ax.set_xticks([])
        ax.set_title(sid, fontsize=8, loc="left", color="black", pad=3)
        for spine in ax.spines.values():
            spine.set_color("black")

        ax.set_xticks(chrom_mid)
        ax.set_xticklabels(
            [c.replace("chr", "") for c in chrom_order],
            fontsize=7,
            color="black"
        )
        ax.set_xlabel("Chromosomes", fontsize=8, color="black")

        for spine in ax.spines.values():
            spine.set_color("black")

    fig.legend(
        *axes[0].get_legend_handles_labels(), loc="upper right",
        fontsize=8, ncol=3, frameon=False,
    )
    fig.suptitle(
        "Per-sample Sadeh local background rate (genome / chromosome / 2 Mb)",
        fontsize=11, color="black",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"[done] -> {args.out}")


if __name__ == "__main__":
    main()
