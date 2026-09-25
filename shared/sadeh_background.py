"""Sadeh local-background-rate model — the single shared implementation.

Three-level hierarchy (genome -> chromosome -> 2 Mb region), trimmed-mean
`fit_noise`. Used by BOTH:
  - build_tissue_markers_sadeh.py (reference side): background-subtracts each
    tissue replicate's raw counts before marker selection.
  - NNLS_deconvolution/scripts/deconvolve.py (patient side): background-
    subtracts each test sample's raw counts before the NNLS solve, so y is
    computed the same way A was.

Both consumers must share this ONE implementation rather than independently
ported copies — otherwise a change to the model (e.g. adding GC-bias
correction) could land in one copy and not the other, silently creating a
scale/bias mismatch between the reference matrix A and the sample vector y.

Kept dependency-light (numpy + pandas only) and free of any snakemake globals
so it can be imported directly by either project.
"""
from __future__ import annotations

import io
import subprocess

import numpy as np
import pandas as pd

# ── Sadeh constants ──────────────────────────────────────────────────────────
_WS_REGION = 2_000_000   # bp — neighbourhood for local rate fitting
_WS_BG_TRIM = 0.95       # Sadeh fitNoise percentile trim
_WS_MIN_BG = 20          # min windows to fit a regional/chromosome rate

# ── GC-bias correction constants ─────────────────────────────────────────────
_GC_N_BINS = 10          # deciles


def region_key_2mb(chrom: np.ndarray, start: np.ndarray) -> np.ndarray:
    """Build the ``chr:<start//2Mb>`` region label used by the regional tier."""
    return np.char.add(
        np.char.add(np.asarray(chrom, dtype=str), ":"),
        (np.asarray(start, dtype=np.int64) // _WS_REGION).astype(str),
    )


def fit_noise(densities: np.ndarray) -> float:
    """Trimmed-mean background rate over finite, non-negative densities."""
    x = densities[np.isfinite(densities) & (densities >= 0)]
    if len(x) < _WS_MIN_BG:
        return float("nan")
    x = x[x <= np.quantile(x, _WS_BG_TRIM)]
    return float(np.mean(x)) if len(x) > 0 else float("nan")


def build_bg_model(bg_counts, bg_len_kb, bg_chr, bg_region2mb, hierarchy) -> dict:
    """Fit genome / per-chromosome / per-2Mb-region background rates (density
    per kb) from the sample's autosomal BACKGROUND-tile counts."""
    density = bg_counts / bg_len_kb
    rate_g = fit_noise(density)
    rate_chr, rate_reg = {}, {}
    if hierarchy in ("chromosome", "regional"):
        for ch in np.unique(bg_chr):
            r = fit_noise(density[bg_chr == ch])
            rate_chr[ch] = r if np.isfinite(r) else rate_g
    if hierarchy == "regional":
        for rk in np.unique(bg_region2mb):
            r = fit_noise(density[bg_region2mb == rk])
            ch = rk.rsplit(":", 1)[0]
            fb = rate_chr.get(ch, rate_g)
            rate_reg[rk] = r if np.isfinite(r) else fb
    return {"genome": rate_g, "chr": rate_chr, "region": rate_reg}


def expected_density_vec(region_key, chr_key, model, hierarchy) -> np.ndarray:
    """Look up the local background density (per kb) for each signal window,
    falling regional -> chromosome -> genome where a finer tier is missing."""
    n = len(region_key)
    out = np.full(n, np.nan)
    if hierarchy == "regional":
        out = pd.Series(region_key).map(model["region"]).to_numpy(dtype=float)
    if hierarchy in ("regional", "chromosome"):
        missing = np.isnan(out)
        if missing.any():
            chr_fill = pd.Series(chr_key).map(model["chr"]).to_numpy(dtype=float)
            out = np.where(missing, chr_fill, out)
    missing = np.isnan(out)
    if missing.any():
        out = np.where(missing, model["genome"], out)
    out = np.where(np.isfinite(out), out, 1e-6)
    return np.maximum(out, 1e-6)


# ── GC-bias correction ───────────────────────────────────────────────────────
# Two separable bias sources, composed rather than conflated:
#   spatial hierarchy (above)  = chromatin-domain effect: is this neighbourhood
#     generally open or closed. Unchanged by anything below.
#   GC correction (here)       = library-prep/sequencing effect: PCR/capture
#     efficiency varies with fragment GC content, roughly constant within one
#     sequencing library, independent of genomic neighbourhood.
# Fit on the RESIDUAL (observed_density / already-fitted spatial rate), not on
# raw density vs GC directly -- fitting on raw density would re-explain some of
# the variance the spatial hierarchy already accounts for (double counting).
# Fit per replicate / per sample, never pooled across libraries: GC bias is a
# property of that specific sequencing run, not a universal constant.
def compute_gc_content(bed_path: str, fasta_path: str) -> "pd.Series":
    """GC fraction (0-1) per region in `bed_path` (region_id must be the 4th
    column), via `bedtools nuc`. Returns a pandas Series (region_id index),
    not a dict -- every caller already does `series.map(gc_map)` /
    `len(gc_map)`, both of which work identically on a Series, but this
    avoids building a multi-million-entry Python dict via a per-line Python
    loop (real memory overhead per entry, and a GIL-holding loop) in favor of
    one vectorized C-level parse."""
    out = subprocess.run(
        ["bedtools", "nuc", "-fi", fasta_path, "-bed", bed_path],
        capture_output=True, text=True, check=True,
    ).stdout
    df = pd.read_csv(io.StringIO(out), sep="\t", skiprows=1, header=None,
                     usecols=[3, 5], names=["region_id", "gc"])
    return df.set_index("region_id")["gc"]


def fit_gc_correction(residual: np.ndarray, gc_frac: np.ndarray,
                       n_bins: int = _GC_N_BINS) -> dict:
    """Bin background-tile residuals (observed_density / spatial-expected
    density) by GC-fraction decile; trimmed-mean per bin (same style as
    `fit_noise`). Returns a model consumed by `gc_correction_vec`.

    A well-behaved residual centers near 1.0 (the spatial rate is itself a
    trimmed mean of these same tiles); the GC-binned deviation from 1.0 is the
    actual bias signal. Falls back to a flat factor of `fallback` (itself
    falling back to 1.0 = no correction) wherever a bin has too few tiles to
    fit reliably, or when there isn't enough background data overall to
    support `n_bins` deciles.
    """
    finite = np.isfinite(residual) & np.isfinite(gc_frac) & (residual >= 0)
    residual, gc_frac = residual[finite], gc_frac[finite]

    fallback = fit_noise(residual) if len(residual) >= _WS_MIN_BG else float("nan")
    fallback = fallback if np.isfinite(fallback) and fallback > 0 else 1.0

    if len(residual) < _WS_MIN_BG * n_bins:
        return {"edges": np.array([0.0, 1.0]), "values": np.array([fallback]),
                "fallback": fallback}

    edges = np.quantile(gc_frac, np.linspace(0.0, 1.0, n_bins + 1))
    edges[0], edges[-1] = 0.0, 1.0
    edges = np.unique(edges)                    # guard against degenerate ties
    bin_idx = np.clip(np.digitize(gc_frac, edges[1:-1], right=True), 0, len(edges) - 2)

    values = np.full(len(edges) - 1, fallback)
    for b in range(len(edges) - 1):
        in_bin = bin_idx == b
        v = fit_noise(residual[in_bin]) if in_bin.sum() >= _WS_MIN_BG else float("nan")
        values[b] = v if np.isfinite(v) and v > 0 else fallback

    return {"edges": edges, "values": values, "fallback": fallback}


def gc_correction_vec(gc_frac: np.ndarray, gc_model: dict) -> np.ndarray:
    """Vectorised lookup of the GC correction factor for each window's GC
    fraction; values outside the fitted bin range clamp to the nearest edge."""
    edges, values = gc_model["edges"], gc_model["values"]
    idx = np.clip(np.digitize(gc_frac, edges[1:-1], right=True), 0, len(values) - 1)
    out = values[idx]
    return np.where(np.isfinite(out) & (out > 0), out, gc_model["fallback"])
