#!/usr/bin/env python3
"""
build_tissue_markers_sadeh.py  —  H3K27ac tissue-specific marker BED files,
Sadeh local-background-subtracted variant of build_tissue_markers.py.

Pipeline
--------
1.  Read Windows.csv; keep rows with annotation type in --window-types
    (default TSS, EXTRA_TSS, Enhancer). Collapse duplicate coordinates to the
    highest-priority type, exactly as build_tissue_markers.py does.
    Separately load BACKGROUND-type rows on autosomes (chr1-22) as the
    per-sample background-rate tile pool (Sadeh's own window universe already
    partitions the genome into TSS/Enhancer/Flanking/Background — no need for
    a separately curated off-target BED).
2.  For each tissue replicate (raw tagAlign.gz or BAM — see TISSUE_REPLICATES
    below), count fragment midpoints once over a combined BED of signal
    windows + background tiles (bedtools intersect -c), fit a per-replicate
    3-level Sadeh background-rate model (genome -> chr -> 2 Mb; ported from
    MS-cfEpiMap/scripts/module4_deconvolution.py), and subtract:
        S(w, replicate) = max(C(w) - local_rate(w) * width_kb(w), 0)
3.  Average S(w, replicate) across a tissue's replicates -> S(w, tissue).
    This replaces bedgraph-averaging of pre-merged tissues: each replicate
    gets its own background fit before averaging, which is the statistically
    valid way to pool count data across biological replicates.
4.  Quantile-normalise S across tissues, restricted to windows overlapping
    housekeeping_seed.bed (the Sadeh-paper "housekeeping genes" role — NOT
    constitutive_anchors.bed, which is derived from the same reference peaks
    being normalised here and would be circular). Fit a per-tissue scale
    factor via robust regression (reused from build_tissue_markers.py) and
    apply it to the full matrix.
5.  Threshold filtering + top-N SNR selection (reused from
    build_tissue_markers.py). Defaults are permissive (--high 0 --low inf)
    since the old thresholds were tuned for raw-pileup units that no longer
    apply here — --top-n is the primary selector until re-tuned.
6.  Write one gzipped BED per tissue + summary TSV + SNR-curve plot, same
    layout as build_tissue_markers.py.

Notes
-----
- Roadmap tagAlign files are single-end 36 bp tags; Breast/Prostate/Neutrophil
  are BAM (also single-end — `bwa samse` per the @PG header), converted with
  `bedtools bamtobed` directly (no samtools dependency needed).
- Prostate's two donors each have their raw reads split across two BAM
  part-files (different ENCFF read-accession pools) — those are concatenated
  into ONE replicate per donor, not treated as two replicates.
- CNS single-cell tracks (Astrocyte/Microglia/Oligodendrocyte/Neuron/OPC/VLMC)
  are intentionally out of scope — only bulk Brain (7 region replicates) is
  included, matching what's already active in build_tissue_markers.py today.
- Neutrophil has no confirmed local raw hg19 alignments yet; fill in
  NEUTROPHIL_BAM_FILES below once the correct ENCODE accession is confirmed.
  Until then the tissue is skipped with a warning.
- Fragment-midpoint extraction and per-region counts are cached under
  --cache-dir so re-runs (e.g. to retune --top-n) don't re-derive them.

Usage
-----
    python build_tissue_markers_sadeh.py [--cpus N] [--high THR] [--low THR]
        [--target-fdr Q] [--n-permutations N] [--top-n N]
        [--window-types TSS,EXTRA_TSS,Enhancer]
        [--housekeeping-bed PATH] [--bg-hierarchy regional|chromosome|global]
        [--cache-dir DIR] [--out DIR]
"""

import argparse
import gzip
import io
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Reuse quantile normalisation / robust-regression helpers instead of
# duplicating them — they're generic (operate on a windows x tissues matrix)
# and unrelated to the raw-vs-background-subtracted distinction.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_tissue_markers import quantile_normalize, fit_rlm_slope  # noqa: E402
# Sadeh background model is shared verbatim with NNLS_deconvolution/scripts/
# deconvolve.py (see shared/sadeh_background.py's docstring for why: A and y
# must be background-subtracted by the exact same code path).
from shared.sadeh_background import (  # noqa: E402
    fit_noise, build_bg_model, expected_density_vec,
    compute_gc_content, fit_gc_correction, gc_correction_vec,
)

# ── Paths ───────────────────────────────────────────────────────────────────
BASE          = Path("/cfs/klemming/projects/supr/uppstore2017150/anita")
CFCHROM_DIR   = BASE / "cfChromatin_in_MS"
WINDOWS_CSV   = BASE / "Sadeh_markers/Windows.csv"
TAGALIGN_DIR  = BASE / "epiroadmap/raw"
BAM_DIR       = BASE / "epiroadmap/bam"
HOUSEKEEPING_BED_DEFAULT = (
    CFCHROM_DIR / "MS-cfEpiMap/reference/rre/housekeeping_seed.bed"
)
FR_DIR = BASE / "Human_Brain_Epi_Atlas/cluster_fragments"
# Blood-open veto (see build_blood_veto_mask). PBMC broad peaks used as-is;
# neutrophil narrow peaks are slopped to approximate broad coverage — both
# confirmed hg19. Neutrophils are granulocytes, absent from PBMC (Ficoll), so
# both compartments are needed to cover the dominant cfDNA sources.
PBMC_BROAD_BED_DEFAULT   = BASE / "epiroadmap/broad_peaks/PBMC_broad_peaks.bed.gz"
MEGAKARYOCYTE_BROAD_BED_DEFAULT =BASE / "epiroadmap/broad_peaks/Megakaryocyte_MEP_broad_peaks.bed.gz"
NEUTROPHIL_PEAKS_DEFAULT = (
    CFCHROM_DIR / "MS-cfEpiMap/reference/rre/raw/Neutrophil_peaks.bed.gz"
)
CHROM_SIZES_DEFAULT      = BASE / "hg19/hg19.chrom.sizes"
FASTA_DEFAULT            = BASE / "hg19/ucsc_hg19.fa"

# ── Annotation type priority (higher = preferred when deduplicating) ────────
ANN_PRIORITY = {"Enhancer": 0, "EXTRA_TSS": 1, "TSS": 2}
SELECTABLE_TYPES = set(ANN_PRIORITY)

AUTOSOME_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2])$")

# ── Tissue -> replicate raw-file map ─────────────────────────────────────────
# tagAlign-backed tissues: name -> list of Roadmap EIDs (one EID = one
# replicate). Filenames follow raw/{EID}-H3K27ac.tagAlign.gz, same mapping
# used by create_roadmap_bedgraphs.py / create_brain_bedgraph.py.
ROADMAP_EIDS = {
    "adipose":           ["E063"],
    "b_cell":            ["E032"],
    "epithelial":        ["E055", "E056", "E058", "E061"],
    "liver":             ["E066"],
    #"lung":              ["E128"],
    "ovary":             ["E097"],
    "nk":                ["E046"],
    "monocyte":          ["E029"],
    #"skin":              ["E127"],
    #"placenta":          ["E091"],
    "gi_smooth_muscle":  ["E111"],
    "skeletal_muscle_f": ["E108"],
    "t_cell":            ["E041", "E042", "E044", "E038", "E037", "E047", "E048"],
    "brain":             ["E067", "E068", "E069", "E071", "E072", "E073", "E074"],
    "digestive":         ["E106", "E102", "E094"],
    "heart":             ["E104", "E105", "E065"],
    #"pancreas":          ["E098", "E087"],
    "blood_vessel":       ["E122"]
}

# BAM-backed tissues: name -> list of (replicate_id, [bam filenames]).
# A replicate with >1 filename means the donor's reads are split across BAM
# part-files that must be concatenated (NOT separate replicates) — this is
# the case for both Prostate donors.
BREAST_BAM_FILES = [
    ("Breast_51", ["Breast_51.bam"]),
    ("Breast_53", ["Breast_53.bam"]),
]
PROSTATE_BAM_FILES = [
    ("Prostate_37", ["Prostate_37_1.bam", "Prostate_37_2.bam"]),
    ("Prostate_54", ["Prostate_54_1.bam", "Prostate_54_2.bam"]),
]
# Neutrophil H3K27ac (ENCODE ENCSR267YXV), both BAMs confirmed hg19.
# Neutrophil_1.bam and Neutrophil_2.bam each independently pool their OWN
# distinct set of 4 source FASTQs (different ENCFF accessions per BAM) -> these
# are the experiment's two BIOLOGICAL replicates, not part-files of one sample.
# They must be fit + normalized separately then averaged (one tuple each), NOT
# concatenated into a single ~2x-depth pseudo-replicate (which inflated
# neutrophil's marker count via depth imbalance in the cross-tissue contest).
NEUTROPHIL_BAM_FILES = [
    ("Neutrophil_1", ["Neutrophil_1.bam"]),
    ("Neutrophil_2", ["Neutrophil_2.bam"]),
]

# Liver: the two BAMs come from disjoint source FASTQs (ENCFF778YPM vs
# ENCFF667AVG+ENCFF343XGX) at very different depths (41M vs 82M) -> two
# biological replicates. Kept separate so each is fit + normalized on its own
# depth and the S(w) columns are averaged; concatenating gave a ~122M-read
# pseudo-replicate (deepest track in the panel) that won spurious markers.
#LIVER_BAM_FILES = [
    #("Liver_1", ["ENCFF998QZR.bam"]),
    #("Liver_2", ["ENCFF774ZCQ.bam"]),
#]

# Sk_muscle: two replicates (disjoint FASTQs ENCFF000BLI vs ENCFF000BLJ). Kept
# separate for consistency, but note both are shallow (6.5M / 9.6M reads) so
# per-replicate Sadeh background fits lean on the chr/global hierarchy fallback;
# if the fit diagnostic degrades, reverting to a single pooled ~16M replicate is
# the reasonable alternative here since pooling does not inflate this tissue.
#SK_MUSCLE_BAM_FILES = [
    #("Sk_muscle_1", ["ENCFF850HGQ.bam"]),
    #("Sk_muscle_2", ["ENCFF273XJE.bam"]),
#]

LUNG_BAM_FILES = [("Lung", ["Lung.bam"])]

MEP_BAM_FILES = [("Megakaryocyte_MEP", ["Megakaryocyte_MEP.bam"])]

SKIN_BAM_FILES = [("Skin", ["Skin_ENCFF753TDL.bam"])]

BAM_TISSUES = {
    "breast":     BREAST_BAM_FILES,
    "prostate":   PROSTATE_BAM_FILES,
    "neutrophil": NEUTROPHIL_BAM_FILES,
    #"liver": LIVER_BAM_FILES,
    #"skeletal_muscle_f": SK_MUSCLE_BAM_FILES,
    "megakaryocyte_mep": MEP_BAM_FILES,
    "lung": LUNG_BAM_FILES,
    "skin": SKIN_BAM_FILES
}

#ASTROCYTE_FR_FILES = [("Astrocyte", ["Astrocyte_fragments.tsv.gz"])]
#MICROGLIA_FR_FILES = [("Microglia", ["Microglia_fragments.tsv.gz"])]
#NEURON_FR_FILES = [("Neuron", ["Neuron_fragments.tsv.gz"])]
#OLIGODENDROCYTE_FR_FILES = [("Oligodendrocyte", ["Oligodendrocyte_fragments.tsv.gz"])]
#OPC_FR_FILES = [("OPC", ["OPC_fragments.tsv.gz"])]
BRAIN_PB_FILES = [("brain_pb", ["merged_fragments.tsv.gz"])]

FRAGMENT_TISSUES = {
    #"brain_pb" : BRAIN_PB_FILES,
    #"astrocyte": ASTROCYTE_FR_FILES,
    #"microglia": MICROGLIA_FR_FILES,
    #"neuron": NEURON_FR_FILES,
    #"oligodendrocyte": OLIGODENDROCYTE_FR_FILES,
    #"opc": OPC_FR_FILES
}

# Hematopoietic tissues exempt from the blood-open veto — a blood-open window is
# a legitimate marker for these. All OTHER tissues get blood-overlapping
# candidate windows dropped before marker selection. Note blood_vessel is NOT
# here: endothelium is a real tissue-of-origin and should be vetoed vs. blood.
BLOOD_TISSUES = {"neutrophil", "monocyte", "b_cell", "t_cell", "nk", "megakaryocyte_mep"}

# ── Sadeh background-model constants (ported from module4_deconvolution.py) ─
_WS_REGION   = 2_000_000    # bp — neighbourhood for local rate fitting
_WS_BG_TRIM  = 0.95         # Sadeh fitNoise percentile trim
_WS_MIN_BG   = 20           # min windows to fit a regional/chromosome rate


def run(cmd: str) -> str:
    result = subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True,
        executable="/bin/bash",
    )
    return result.stdout


# ── Replicate resolution ─────────────────────────────────────────────────────

def build_tissue_replicates() -> dict[str, list[dict]]:
    """tissue -> [{'id':..., 'filetype': 'tagalign'|'bam', 'paths': [Path,...]}]"""
    replicates: dict[str, list[dict]] = {}

    for tissue, eids in ROADMAP_EIDS.items():
        reps = []
        for eid in eids:
            p = TAGALIGN_DIR / f"{eid}-H3K27ac.tagAlign.gz"
            if not p.exists():
                print(f"  SKIP replicate {eid} ({tissue}): {p.name} not found.")
                continue
            reps.append({"id": eid, "filetype": "tagalign", "paths": [p]})
        if reps:
            replicates[tissue] = reps

    for tissue, members in BAM_TISSUES.items():
        if not members:
            print(f"  SKIP {tissue}: no BAM files configured yet.")
            continue
        reps = []
        for rep_id, filenames in members:
            paths = [BAM_DIR / f for f in filenames]
            missing = [p for p in paths if not p.exists()]
            if missing:
                print(f"  SKIP replicate {rep_id} ({tissue}): "
                      f"missing {[p.name for p in missing]}")
                continue
            reps.append({"id": rep_id, "filetype": "bam", "paths": paths})
        if reps:
            replicates[tissue] = reps

    for tissue, members in FRAGMENT_TISSUES.items():
        if not members:
            print(f"  SKIP {tissue}: no fragment files configured yet.")
            continue
        reps = []
        for rep_id, filenames in members:
            paths = [FR_DIR / f for f in filenames]
            missing = [p for p in paths if not p.exists()]
            if missing:
                print(f"  SKIP replicate {rep_id} ({tissue}): "
                      f"missing {[p.name for p in missing]}")
                continue
            reps.append({"id": rep_id, "filetype": "fragments", "paths": paths})
        if reps:
            replicates[tissue] = reps

    return replicates


# ── Window loading (Windows.csv) ─────────────────────────────────────────────

def load_signal_windows(window_types: list[str]) -> pd.DataFrame:
    """Load + dedup Windows.csv rows matching window_types (highest-priority
    annotation wins on coordinate collisions), mirroring build_tissue_markers.py
    Step 1."""
    windows_raw = pd.read_csv(
        WINDOWS_CSV, header=None,
        names=["chr", "start", "end", "gene", "tissue_ann", "type"],
        dtype={"chr": str, "start": int, "end": int,
               "gene": str, "tissue_ann": str, "type": str},
    )
    keep = windows_raw["type"].isin(window_types)
    windows_filt = windows_raw[keep].copy()

    windows_filt["_pri"] = windows_filt["type"].map(ANN_PRIORITY)
    windows_filt = (
        windows_filt
        .sort_values("_pri", ascending=False)
        .drop_duplicates(subset=["chr", "start", "end"], keep="first")
        .drop(columns="_pri")
        .sort_values(["chr", "start"])
        .reset_index(drop=True)
    )
    windows_filt["region_id"] = (
        "sig|" + windows_filt["chr"] + ":"
        + windows_filt["start"].astype(str) + "-"
        + windows_filt["end"].astype(str)
    )
    windows_filt["length_kb"] = (
        (windows_filt["end"] - windows_filt["start"]) / 1000.0
    )
    windows_filt["region_2mb"] = (
        windows_filt["chr"] + ":"
        + (windows_filt["start"] // _WS_REGION).astype(str)
    )
    # category dtype after the string-concat steps above (which need plain
    # strings) -- chr/type each have a handful of distinct values across
    # 4M+ rows, so this cuts real memory with no behavior change downstream
    # (.to_numpy()/comparisons/.map() all work the same on category columns).
    windows_filt["chr"] = windows_filt["chr"].astype("category")
    windows_filt["type"] = windows_filt["type"].astype("category")
    return windows_filt


def load_background_tiles() -> pd.DataFrame:
    """Load Windows.csv BACKGROUND rows on autosomes as the background-rate
    tile pool (replaces module4's curated off-target BED)."""
    windows_raw = pd.read_csv(
        WINDOWS_CSV, header=None,
        names=["chr", "start", "end", "gene", "tissue_ann", "type"],
        dtype={"chr": str, "start": int, "end": int,
               "gene": str, "tissue_ann": str, "type": str},
    )
    bg = windows_raw[windows_raw["type"] == "BACKGROUND"].copy()
    bg = bg[bg["chr"].apply(lambda c: bool(AUTOSOME_RE.match(c)))]
    bg = bg.sort_values(["chr", "start"]).reset_index(drop=True)
    bg["region_id"] = (
        "bg|" + bg["chr"] + ":" + bg["start"].astype(str) + "-"
        + bg["end"].astype(str)
    )
    bg["length_kb"] = (bg["end"] - bg["start"]) / 1000.0
    bg["region_2mb"] = bg["chr"] + ":" + (bg["start"] // _WS_REGION).astype(str)
    bg["chr"] = bg["chr"].astype("category")
    bg["type"] = bg["type"].astype("category")
    return bg


# ── Fragment extraction + counting (per replicate, cached) ──────────────────

def get_midpoint_bed(rep: dict, cache_dir: Path) -> Path:
    """Sorted BED of read midpoints for one replicate, cached to disk.
    Independent of window selection, so changing --window-types never forces
    re-extraction. filetype: 'tagalign' | 'bam' | 'fragments'."""
    out_gz = cache_dir / f"{rep['id']}.midpoints.bed.gz"
    if out_gz.exists():
        return out_gz

    if rep["filetype"] == "tagalign":
        (path,) = rep["paths"]
        raw_cmd = f"zcat {path}"
        midpoint_awk = (
            "awk 'BEGIN{OFS=\"\\t\"}{mid=int(($2+$3)/2); print $1,mid,mid+1}'"
        )
    elif rep["filetype"] == "bam":
        raw_cmd = "(" + " && ".join(
            f"bedtools bamtobed -i {p}" for p in rep["paths"]
        ) + ")"
        midpoint_awk = (
            "awk 'BEGIN{OFS=\"\\t\"}{mid=int(($2+$3)/2); print $1,mid,mid+1}'"
        )
    elif rep["filetype"] == "fragments":
        # 10x Cell Ranger ATAC-style fragments.tsv.gz, pre-aggregated per cell
        # type: chrom, start, end, barcode, count (col 5 = number of PCR-
        # duplicate reads supporting this exact fragment). Two things this
        # format needs that tagalign/bam don't:
        #   - a '#'-prefixed metadata header block to skip;
        #   - the count column EXPANDED into that many midpoint copies, not
        #     collapsed to 1 -- tagAlign/BAM lines are one-row-per-READ
        #     (duplicates included), so to keep read-depth comparable across
        #     tissues (the Sadeh background rate + cross-tissue brightness
        #     scaling both assume "count" means read depth) a fragment with
        #     count=5 must contribute 5 midpoints, matching what 5 separate
        #     tagAlign/BAM lines at that position would have contributed. The
        #     barcode column (4th field) is unused: these files are already
        #     pooled per cell type, no per-cell grouping needed here.
        (path,) = rep["paths"]
        raw_cmd = f"zcat {path}"
        midpoint_awk = (
            "awk 'BEGIN{OFS=\"\\t\"} /^#/{next} "
            "{mid=int(($2+$3)/2); for(i=0;i<$5;i++) print $1,mid,mid+1}'"
        )
    else:
        raise ValueError(f"Unknown filetype: {rep['filetype']}")

    tmp_out = out_gz.with_suffix(".tmp.gz")
    # pipefail so a failure anywhere upstream (corrupt/empty input) raises
    # here instead of being swallowed by gzip's own exit code.
    # Canonical-chromosome filter: fragment files from 10x hg38 references
    # keep alt/unplaced scaffolds under bare GenBank accessions (e.g.
    # GL000009.2, no "chr" prefix) instead of UCSC's chr*_alt/_random
    # naming. combined_bed never has these accessions, and the mismatched
    # naming convention breaks `bedtools intersect -sorted`'s chromosome
    # ordering (hard failure, not just a missed match), so drop them here.
    run(
        f"set -o pipefail; {raw_cmd} | {midpoint_awk} | "
        f"awk '$1 ~ /^chr([0-9]+|X|Y)$/' | "
        f"sort -k1,1 -k2,2n | gzip > {tmp_out}"
    )
    n_lines = int(run(f"zcat {tmp_out} | wc -l").strip() or 0)
    if n_lines == 0:
        tmp_out.unlink()
        raise RuntimeError(
            f"[{rep['id']}] extracted 0 reads from {rep['paths']} — "
            "check the source file(s) are valid, non-empty alignments.")
    tmp_out.rename(out_gz)
    return out_gz


def count_replicate(rep_id: str, midpoint_bed: Path, combined_bed: Path,
                     wt_sig: str, cache_dir: Path) -> pd.Series:
    """region_id -> count for one replicate, via bedtools intersect -c
    against the combined signal+background BED (mirrors module4's
    _count_per_region_sample).

    Returns a pandas Series (region_id index), not a dict. A ~4.6M-entry
    Python dict plus a per-element `.get()` loop to rebuild bg/sig arrays
    (in process_replicate) costs real memory (dict overhead per entry) and
    CPU (a pure-Python loop, done for EVERY replicate, including cache
    hits -- the cache saves the bedtools call, not this reconstruction).
    It also holds the GIL for the whole loop, which quietly serializes the
    ThreadPoolExecutor workers in Step 4 during exactly this step. A Series
    supports the same downstream lookup via a single vectorized
    `.reindex(...)` call instead (see process_replicate) -- same result,
    computed in C, and releases the GIL while doing it.
    """
    cache_path = cache_dir / f"{rep_id}.{wt_sig}.counts.tsv.gz"
    if cache_path.exists():
        df = pd.read_csv(cache_path, sep="\t", names=["region_id", "count"])
        return df.set_index("region_id")["count"]

    # -sorted streams both inputs chromosome-by-chromosome instead of loading
    # the entire (tens-of-millions-of-reads) -b file into an in-memory
    # interval tree — both inputs are already sort -k1,1 -k2,2n, so this is
    # safe and cuts peak memory drastically under parallel workers.
    out = run(
        f"bedtools intersect -a {combined_bed} -b {midpoint_bed} -c -sorted"
    )
    # bedtools' own columns are chr,start,end,region_id,...,count(last) --
    # combined_bed is exactly 4 columns so count is column index 4.
    df = pd.read_csv(io.StringIO(out), sep="\t", header=None,
                     usecols=[3, 4], names=["region_id", "count"])

    tmp_cache = cache_path.with_suffix(".tmp.gz")
    df.to_csv(tmp_cache, sep="\t", header=False, index=False, compression="gzip")
    tmp_cache.rename(cache_path)
    return df.set_index("region_id")["count"]


# Sadeh background model (fit_noise, build_bg_model, expected_density_vec) now
# lives in shared/sadeh_background.py — imported above, not duplicated here.


# ── Per-replicate worker ─────────────────────────────────────────────────────

def process_replicate(rep: dict, tissue: str, combined_bed: Path, wt_sig: str,
                       sig_df: pd.DataFrame, bg_arrays: dict, hierarchy: str,
                       cache_dir: Path, gc_correction: bool = True
                       ) -> tuple[str, str, np.ndarray]:
    midpoint_bed = get_midpoint_bed(rep, cache_dir)
    counts = count_replicate(rep["id"], midpoint_bed, combined_bed, wt_sig,
                              cache_dir)

    bg_counts = counts.reindex(bg_arrays["region_id"], fill_value=0).to_numpy(dtype=float)
    model = build_bg_model(bg_counts, bg_arrays["length_kb"],
                            bg_arrays["chr"], bg_arrays["region_2mb"],
                            hierarchy)

    sig_counts = counts.reindex(sig_df["region_id"].to_numpy(), fill_value=0).to_numpy(dtype=float)
    exp_density = expected_density_vec(
        sig_df["region_2mb"].to_numpy(), sig_df["chr"].to_numpy(),
        model, hierarchy)

    # GC-bias correction (residual, fit on THIS replicate's own bg tiles — see
    # shared/sadeh_background.py's docstring for why it's a residual, not a
    # replacement, for the spatial hierarchy above).
    if gc_correction:
        bg_len_kb = bg_arrays["length_kb"]
        exp_bg = expected_density_vec(
            bg_arrays["region_2mb"], bg_arrays["chr"], model, hierarchy)
        bg_residual = (bg_counts / bg_len_kb) / np.maximum(exp_bg, 1e-9)
        gc_model = fit_gc_correction(bg_residual, bg_arrays["gc"])
        gc_factor = gc_correction_vec(sig_df["gc"].to_numpy(), gc_model)
    else:
        gc_factor = np.ones(len(sig_df))

    exp_count = exp_density * gc_factor * np.maximum(sig_df["length_kb"].to_numpy(), 1e-3)

    signal = np.maximum(sig_counts - exp_count, 0.0).astype(np.float32)
    print(f"    done: {tissue}/{rep['id']}  "
          f"(genome rate={model['genome']:.4g} reads/kb)")
    return tissue, rep["id"], signal, sig_counts.astype(np.float32)


def mean_variance_plot(count_matrix: np.ndarray, out_path: Path, title: str,
                       unit_label: str) -> None:
    """Per-window mean-vs-variance scatter across columns (samples/replicates).

    Diagnoses count overdispersion. Poisson predicts var == mean (red y=x);
    a negative-binomial model predicts var = mean + alpha*mean^2. A single
    global alpha is fit by regressing (var - mean) on mean^2 through the origin
    and drawn as the trend — the same overdispersion the per-replicate Sadeh
    model summarises as alpha_s, but here across the panel rather than within
    background tiles of one replicate.
    """
    M = np.asarray(count_matrix, dtype=float)
    if M.ndim != 2 or M.shape[1] < 2:
        print(f"  [mean-variance] need >=2 columns, got shape {M.shape}; skipping.")
        return
    m = M.mean(axis=1)
    v = M.var(axis=1, ddof=1)
    keep = m > 0
    m, v = m[keep], v[keep]
    if m.size == 0:
        print("  [mean-variance] no non-zero windows; skipping.")
        return
    # NB overdispersion: var = mean + alpha*mean^2  ->  regress (v-m) on m^2 thru 0.
    m2 = m ** 2
    denom = float(np.sum(m2 ** 2))
    alpha = max(float(np.sum((v - m) * m2) / denom), 0.0) if denom > 0 else 0.0

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(m, v, s=6, alpha=0.25, edgecolors="none", color="#3b6fb0",
               label=f"windows (n={m.size:,})")
    grid = np.geomspace(max(float(m.min()), 1e-3), float(m.max()), 200)
    ax.plot(grid, grid, "r--", lw=1.5, label="Poisson (var = mean)")
    ax.plot(grid, grid + alpha * grid ** 2, "k-", lw=1.5,
            label=f"NB trend (α = {alpha:.3g})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"Per-window mean count  ({unit_label})")
    ax.set_ylabel("Per-window variance")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  mean-variance plot → {out_path}")


def correlation_heatmap(matrix: np.ndarray, tissues: list[str],
                        out_path: Path, title: str) -> None:
    """Annotated Pearson-correlation heatmap between tissue signature columns.

    `matrix` is windows x tissues (restricted to the marker/universe rows so the
    correlations reflect the actual columns that enter the deconvolution matrix
    A). High off-diagonal values flag collinear tissues that NNLS cannot cleanly
    separate — the empirical basis for deciding whether to lump columns (e.g.
    the lymphocyte question) or drop a redundant tissue.
    """
    M = np.asarray(matrix, dtype=float)
    if M.ndim != 2 or M.shape[1] < 2 or M.shape[0] < 2:
        print(f"  [correlation] need >=2 rows and >=2 columns, got {M.shape}; skipping.")
        return
    corr = np.corrcoef(M, rowvar=False)          # tissues x tissues
    corr = np.nan_to_num(corr, nan=0.0)          # constant columns -> 0 correlation
    n = len(tissues)

    fig, ax = plt.subplots(figsize=(max(7, 0.7 * n), max(6, 0.7 * n)))
    im = ax.imshow(corr, vmin=-1.0, vmax=1.0, cmap="RdBu_r", aspect="equal")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(tissues, rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=8)
    ax.set_yticklabels(tissues, fontsize=8)

    # Number inside every cell; white text on saturated cells for legibility.
    for i in range(n):
        for j in range(n):
            r = corr[i, j]
            ax.text(j, i, f"{r:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(r) > 0.5 else "black")

    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  correlation heatmap → {out_path}")


# ── GC-correction diagnostics (prove the correction is warranted) ────────────

def gc_correction_diagnostics(flat_reps, combined_bed, wt_sig, sig_df,
                              bg_arrays, hierarchy, cache_dir, out_dir,
                              n_bins: int = 10, seed: int = 0) -> None:
    """Evidence that the per-replicate GC-bias correction is both warranted and
    non-harmful — the reference-side analogue of module4's within-sample GC
    checks (residual model + permuted-signature null) that this script
    otherwise lacked while still APPLYING the correction (process_replicate).

    Two diagnostics, written to <out>/qc/:
      1. Residual-vs-GC, before vs after. For every replicate, bin its
         background-tile residual (observed_density / spatial-expected) by GC
         decile. A pre-correction curve that bows away from 1.0 at the GC
         extremes IS the bias the spatial hierarchy misses; the post-correction
         curve (residual / fitted GC factor) should flatten. The INFORMATIVE
         quantity is the AMPLITUDE of the pre-correction curve — post-flatness
         is partly by construction (the factor is the per-decile trimmed mean),
         so it is a correctness check, not proof on its own.
      2. GC×length-matched fake-signature null. Draw background tiles matched to
         the real signal windows' GC×length and compute their pooled log2FE
         against spatial-only vs spatial×GC expected. Fake signatures carry no
         tissue signal, so a |median log2FE| that moves toward 0 after GC
         correction is direct evidence the correction removes inflation rather
         than inventing it. Signal windows are GC-rich TSS/enhancers, so this
         probes exactly the high-GC regime where over/under-correction bites.
    """
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    bg_gc   = np.asarray(bg_arrays["gc"], dtype=float)
    bg_len  = np.asarray(bg_arrays["length_kb"], dtype=float)
    bg_chr  = bg_arrays["chr"]
    bg_r2mb = bg_arrays["region_2mb"]
    bg_rid  = bg_arrays["region_id"]

    if not np.isfinite(bg_gc).any() or np.allclose(bg_gc, 0.0):
        print("  [gc-diag] background GC content unavailable — skipping.")
        return

    # GC decile bins for the residual curve (shared across replicates).
    # Bin/aggregate over ONLY finite, non-gap tiles:
    #   - np.nanquantile (not np.quantile): a single NaN GC would otherwise
    #     turn every decile edge into NaN and collapse the curve to 2 junk bins.
    #   - drop GC == 0 assembly-gap / N-run tiles: not a real low-GC regime, and
    #     their spatial rate collapses to the floor → inf residuals.
    # This masking is DISPLAY-ONLY; it does not touch the production GC fit.
    gc_ok = np.isfinite(bg_gc) & (bg_gc > 0.0)
    edges = np.nanquantile(bg_gc[gc_ok], np.linspace(0.0, 1.0, n_bins + 1))
    edges[0], edges[-1] = 0.0, 1.0
    edges = np.unique(edges)
    n_dec = len(edges) - 1
    dec_idx = np.clip(np.digitize(bg_gc, edges[1:-1], right=True), 0, n_dec - 1)
    dec_gc  = np.array(
        [np.median(bg_gc[(dec_idx == b) & gc_ok]) if ((dec_idx == b) & gc_ok).any()
         else np.nan for b in range(n_dec)])

    # GC×length-matched fake-signature indices. Matching is on genomic GC/length
    # (replicate-independent) so the SAME fake set is reused for every replicate;
    # only the counts differ. Cap the sampled signal set so the O(n) match loop
    # stays fast on a multi-million-row signal universe.
    rng = np.random.default_rng(seed)
    _sig = sig_df.sample(20000, random_state=seed) if len(sig_df) > 20000 else sig_df
    sig_gc  = _sig["gc"].to_numpy(dtype=float)
    sig_len = _sig["length_kb"].to_numpy(dtype=float)

    def _bin_keys(gc, length, gc_b=5, len_b=5):
        ge = np.quantile(bg_gc,  np.linspace(0, 1, gc_b + 1)[1:-1])
        le = np.quantile(bg_len, np.linspace(0, 1, len_b + 1)[1:-1])
        gi = np.clip(np.searchsorted(ge, gc), 0, gc_b - 1)
        li = np.clip(np.searchsorted(le, length), 0, len_b - 1)
        return gi * len_b + li

    _bgk = _bin_keys(bg_gc, bg_len)
    _sgk = _bin_keys(sig_gc, sig_len)
    _bin_to_idx = {int(k): np.where(_bgk == k)[0] for k in np.unique(_bgk)}
    _all_bg = np.arange(len(bg_gc))
    fake_idx = np.array(
        [rng.choice(_bin_to_idx.get(int(k), _all_bg)) for k in _sgk], dtype=int)

    before_curves, after_curves, fake_rows = [], [], []
    eps = 1.0

    for rep, tissue in flat_reps:
        try:
            mid = get_midpoint_bed(rep, cache_dir)
            counts = count_replicate(rep["id"], mid, combined_bed, wt_sig, cache_dir)
        except Exception as exc:
            print(f"  [gc-diag] {tissue}/{rep['id']}: counts unavailable ({exc})")
            continue
        bg_counts = counts.reindex(bg_rid, fill_value=0).to_numpy(dtype=float)
        model  = build_bg_model(bg_counts, bg_len, bg_chr, bg_r2mb, hierarchy)
        exp_bg = expected_density_vec(bg_r2mb, bg_chr, model, hierarchy)
        resid  = (bg_counts / bg_len) / np.maximum(exp_bg, 1e-9)
        gc_model = fit_gc_correction(resid, bg_gc)
        resid_after = resid / np.maximum(gc_correction_vec(bg_gc, gc_model), 1e-9)

        # Per-decile TRIMMED mean (same 95th-pct trim as the fitter's fit_noise),
        # over non-gap tiles whose spatial rate did not hit the 1e-6 floor —
        # floor tiles give inf residuals that would swamp an untrimmed mean.
        curve_ok = gc_ok & np.isfinite(resid) & (exp_bg > 1e-6)
        before_curves.append(np.array(
            [fit_noise(resid[(dec_idx == b) & curve_ok]) for b in range(n_dec)]))
        after_curves.append(np.array(
            [fit_noise(resid_after[(dec_idx == b) & curve_ok]) for b in range(n_dec)]))

        # Fake-signature null (matched to real signal GC×length).
        f_cnt = bg_counts[fake_idx]
        f_len = bg_len[fake_idx]
        f_exp = exp_bg[fake_idx]
        f_gcf = gc_correction_vec(bg_gc[fake_idx], gc_model)
        obs        = float(f_cnt.sum())
        exp_uncorr = float((f_exp * f_len).sum())
        exp_corr   = float((f_exp * f_gcf * f_len).sum())
        fake_rows.append({
            "tissue": tissue, "replicate": rep["id"],
            "fake_log2FE_uncorrected": float(np.log2((obs + eps) / (exp_uncorr + eps))),
            "fake_log2FE_corrected":   float(np.log2((obs + eps) / (exp_corr + eps))),
        })

    if not before_curves:
        print("  [gc-diag] no replicate counts available — skipping.")
        return

    before = np.nanmean(np.vstack(before_curves), axis=0)
    after  = np.nanmean(np.vstack(after_curves), axis=0)
    amp_before = float(np.nanmax(np.abs(before - 1.0)))
    amp_after  = float(np.nanmax(np.abs(after - 1.0)))

    fake_df = pd.DataFrame(fake_rows)
    fake_df.to_csv(qc_dir / "gc_fake_signature_null.tsv", sep="\t", index=False)
    pd.DataFrame({"gc_decile_center": dec_gc,
                  "residual_before": before,
                  "residual_after": after}).to_csv(
        qc_dir / "gc_residual_curve.tsv", sep="\t", index=False)

    med_uncorr = float(fake_df["fake_log2FE_uncorrected"].median())
    med_corr   = float(fake_df["fake_log2FE_corrected"].median())

    # ── Figure 1: residual-vs-GC before/after ────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(dec_gc, before, "o-", color="#c0392b", lw=1.5,
            label=f"before  (max |dev| = {amp_before:.0%})")
    ax.plot(dec_gc, after, "o-", color="#2c7fb8", lw=1.5,
            label=f"after   (max |dev| = {amp_after:.0%})")
    ax.axhline(1.0, ls="--", lw=0.8, color="grey")
    ax.set_xlabel("Background-tile GC fraction (decile center)")
    ax.set_ylabel("Mean residual  (observed density / spatial-expected)")
    ax.set_title("GC-bias correction — background-tile residual\n"
                 "(flat at 1.0 = no residual GC bias)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(qc_dir / "gc_residual_before_after.png", dpi=150)
    plt.close(fig)

    # ── Figure 2: fake-signature null, uncorrected vs corrected ──────────────
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(fake_df["fake_log2FE_uncorrected"],
               fake_df["fake_log2FE_corrected"],
               s=30, alpha=0.8, edgecolors="black", lw=0.4, color="#4878CF")
    _lim = float(np.nanmax(np.abs(fake_df[[
        "fake_log2FE_uncorrected", "fake_log2FE_corrected"]].to_numpy())) or 0.1)
    _lim = max(_lim * 1.1, 0.05)
    ax.plot([-_lim, _lim], [-_lim, _lim], ls="--", lw=0.8, color="grey")
    ax.axhline(0, ls=":", lw=0.8, color="red")
    ax.axvline(0, ls=":", lw=0.8, color="red")
    ax.set_xlim(-_lim, _lim); ax.set_ylim(-_lim, _lim)
    ax.set_xlabel("Fake-signature log2FE — uncorrected")
    ax.set_ylabel("Fake-signature log2FE — GC-corrected")
    ax.set_title("GC×length-matched fake-signature null\n"
                 "(closer to 0 on y than x ⇒ correction removes GC inflation)")
    fig.tight_layout()
    fig.savefig(qc_dir / "gc_fake_signature_null.png", dpi=150)
    plt.close(fig)

    # ── Verdict ──────────────────────────────────────────────────────────────
    print(f"  [gc-diag] residual-vs-GC amplitude: before {amp_before:.0%} → "
          f"after {amp_after:.0%}")
    print(f"  [gc-diag] fake-signature median log2FE: uncorrected {med_uncorr:+.3f} "
          f"→ corrected {med_corr:+.3f}")
    if amp_before >= 0.10 and abs(med_corr) < abs(med_uncorr):
        print("  [gc-diag] VERDICT: GC bias is present and the correction reduces "
              "it → correction is warranted.")
    elif amp_before <= 0.05:
        print("  [gc-diag] VERDICT: GC bias is negligible (<5% at the extremes); "
              "--no-gc-correction is defensible and simpler.")
    else:
        print("  [gc-diag] VERDICT: GC bias is modest; inspect the two QC plots "
              "before committing to the correction.")
    print(f"  [gc-diag] QC → {qc_dir}/gc_residual_before_after.png, "
          f"gc_fake_signature_null.png (+ .tsv)")


# ── Blood-open veto mask ─────────────────────────────────────────────────────

def build_blood_veto_mask(sig_df: pd.DataFrame, broad_beds: list[Path],
                          narrow_beds: list[Path], slop_bp: int,
                          chrom_sizes: Path, veto_frac: float,
                          cache_dir: Path) -> np.ndarray:
    """Boolean array aligned to sig_df rows: True where a signal window is
    >= veto_frac covered by blood-open chromatin.

    The mask = union of broad_beds (used as-is) and narrow_beds (slopped by
    slop_bp/side to approximate broad coverage), then merged so a window spanning
    adjacent peaks is seen as one contiguous overlap (bedtools -f measures
    fraction against a single B feature). Non-blood tissues have their
    blood-overlapping candidate windows dropped before selection, so blood-open
    leakage cannot masquerade as tissue-specific signal in blood-dominated cfDNA.
    Mask is cached (keyed by slop) so threshold-tuning re-runs skip the rebuild;
    change the peak files → clear the cached blood_mask.*.bed.
    """
    if not broad_beds and not narrow_beds:
        return np.zeros(len(sig_df), dtype=bool)

    mask_bed = cache_dir / f"blood_mask.slop{slop_bp}.bed"
    if not mask_bed.exists():
        parts = [f"zcat -f {b} | cut -f1-3" for b in broad_beds]
        parts += [
            f"zcat -f {b} | cut -f1-3 | bedtools slop -b {slop_bp} -g {chrom_sizes}"
            for b in narrow_beds
        ]
        # set -eo pipefail so a failed/empty input aborts instead of silently
        # producing a partial mask that would under-veto.
        run(f"set -eo pipefail; ( {' ; '.join(parts)} ) | "
            f"sort -k1,1 -k2,2n | bedtools merge -i - > {mask_bed}")

    tmp_win = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    for idx, r in enumerate(sig_df.itertuples(index=False)):
        tmp_win.write(f"{r.chr}\t{r.start}\t{r.end}\t{idx}\n")
    tmp_win.close()
    tmp_sorted = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    run(f"sort -k1,1 -k2,2n {tmp_win.name} > {tmp_sorted}")
    Path(tmp_win.name).unlink()

    hit = run(
        f"bedtools intersect -a {tmp_sorted} -b {mask_bed} -u -f {veto_frac}")
    Path(tmp_sorted).unlink()

    veto = np.zeros(len(sig_df), dtype=bool)
    for line in hit.splitlines():
        veto[int(line.split("\t")[3])] = True
    return veto


def snr_denom_floor(matrix: np.ndarray, quantile: float = 0.10) -> float:
    """Robust floor for the SNR denominator (other_max), replacing a fixed
    additive epsilon (1e-6).

    Many windows genuinely have other_max == 0 (a true marker has zero signal
    in unrelated tissues after background subtraction) -- dividing by 1e-6
    there gives a huge but real SNR. The problem is the permutation null:
    independently shuffling ~20 columns across millions of windows makes it
    common for SEVERAL columns to coincidentally land a near-zero value on
    the SAME random window, purely by chance. That coincidence also creates a
    near-zero denominator, and dividing by 1e-6 inflates it into an enormous
    *null* SNR with no biological meaning -- common enough (given millions of
    windows) to make the null's tail just as extreme as real markers',
    swamping every FDR comparison and collapsing the selected panel to ~0
    markers for nearly every tissue (confirmed empirically: n_star_fdr == 0
    for brain, liver, adipose... everything except a handful of exact-zero
    windows with astronomical reported SNR).

    Flooring the denominator at a data-driven "practical detection floor" (a
    low percentile of the matrix's own nonzero values, not an arbitrary tiny
    constant) bounds the ratio sanely in BOTH the observed and permuted
    computations -- applied consistently, so the null-vs-observed comparison
    stays meaningful instead of being dominated by division artifacts.
    """
    nonzero = matrix[matrix > 0]
    if len(nonzero) == 0:
        return 1e-6
    floor = float(np.quantile(nonzero, quantile))
    return floor if floor > 0 else 1e-6


# ── Permutation-null FDR for the variable top-N cutoff ──────────────────────
def permutation_null_snr(norm_matrix: np.ndarray, available: list[str],
                          high: float, low: float,
                          blood_veto_mask: np.ndarray, blood_tissues: set,
                          n_permutations: int, rng: np.random.Generator,
                          snr_floor: float,
                          high_pct: float = 0.0, low_pct: float = 0.0
                          ) -> dict[str, np.ndarray]:
    """Null SNR distribution per tissue under "no true tissue-window
    association": independently permute each tissue COLUMN's values across
    windows (breaks cross-tissue joint structure, preserves each tissue's own
    marginal distribution), then run windows through the EXACT SAME selection
    funnel (focal > high & other_max < low & not blood-vetoed for non-blood
    tissues) used on the real data, so the null is a fair comparison against
    the actual selection procedure -- not just "all windows".

    Cost note: computing "max of all OTHER columns" per tissue naively is
    O(n_tissues^2) per permutation. Instead compute the row-wise top-2 values
    and their argmax ONCE per permutation (O(n_tissues) via np.partition) --
    other_max for tissue j is the row's max unless tissue j WAS that row's
    max, in which case it's the second-highest value.
    """
    tissue_idx = {t: i for i, t in enumerate(available)}
    null_snr: dict[str, list] = {t: [] for t in available}
    n_cols = norm_matrix.shape[1]

    # Per-column quantile thresholds, when enabled. Permuting WITHIN a column
    # leaves that column's marginal untouched, so these are the same for every
    # permutation and for the real matrix -- hoist them out of the loop.
    use_pct_p = high_pct > 0.0 and low_pct > 0.0
    if use_pct_p:
        hi_thr_p = np.array([np.quantile(norm_matrix[:, k], high_pct)
                             for k in range(n_cols)])
        lo_thr_p = np.array([np.quantile(norm_matrix[:, k], low_pct)
                             for k in range(n_cols)])

    for _ in range(n_permutations):
        perm = np.empty_like(norm_matrix)
        for j in range(n_cols):
            perm[:, j] = rng.permutation(norm_matrix[:, j])

        if use_pct_p:
            # WHICH rows sit below each threshold does change per permutation,
            # even though the thresholds themselves do not.
            below_lo_p = perm < lo_thr_p[np.newaxis, :]
            n_above_p  = (~below_lo_p).sum(axis=1)

        argmax_col = np.argmax(perm, axis=1)
        if n_cols >= 2:
            top2 = np.partition(perm, -2, axis=1)[:, -2:]
            top1, second = top2[:, 1], top2[:, 0]
        else:
            top1 = perm[:, 0]
            second = np.zeros_like(top1)

        for tissue in available:
            j = tissue_idx[tissue]
            focal_p = perm[:, j]
            other_max_p = np.where(argmax_col == j, second, top1)

            if use_pct_p:
                # Thresholds are quantiles of each column's own marginal, and
                # the permutation shuffles values WITHIN a column -- so the
                # marginal, and therefore hi_thr/lo_thr, are identical to the
                # unpermuted matrix. Compute them once outside the loop.
                others_ok_p = (n_above_p == 0) | ((n_above_p == 1) & ~below_lo_p[:, j])
                is_marker_p = (focal_p > hi_thr_p[j]) & others_ok_p
            else:
                is_marker_p = (focal_p > high) & (other_max_p < low)
            if tissue not in blood_tissues:
                is_marker_p &= ~blood_veto_mask
            idx_p = np.where(is_marker_p)[0]
            snr_p = focal_p[idx_p] / np.maximum(other_max_p[idx_p], snr_floor)
            null_snr[tissue].append(snr_p)

    return {t: (np.concatenate(v) if v else np.array([])) for t, v in null_snr.items()}


def fdr_controlled_n(observed_snr_desc: np.ndarray, null_snr: np.ndarray,
                      n_permutations: int, target_fdr: float) -> tuple[int, float]:
    """Largest N such that the permutation-estimated FDR at rank N is <=
    target_fdr, where FDR(N) = (expected # null exceedances at the Nth
    observed SNR threshold) / N. FDR is stabilized with a running maximum
    from rank 1 upward (standard step-up guard against noise letting FDR
    fictitiously dip back below target at a later, less-trustworthy rank).
    Returns (N, fdr_at_N) -- (0, 0.0) if no rank clears the target.
    """
    n = len(observed_snr_desc)
    if n == 0:
        return 0, 0.0
    if len(null_snr) == 0:
        # No null data (e.g. permutations disabled) -- can't estimate FDR.
        return n, float("nan")

    null_sorted = np.sort(null_snr)
    n_null_total = len(null_sorted)
    idx = np.searchsorted(null_sorted, observed_snr_desc, side="left")
    null_exceed = (n_null_total - idx) / n_permutations
    ranks = np.arange(1, n + 1)
    fdr = null_exceed / ranks
    fdr_guarded = np.maximum.accumulate(fdr)   # monotone non-decreasing floor

    ok = np.where(fdr_guarded <= target_fdr)[0]
    if len(ok) == 0:
        return 0, float(fdr_guarded[0]) if n > 0 else 0.0
    n_star = int(ok[-1] + 1)
    return n_star, float(fdr_guarded[n_star - 1])


# ── Marker selection (reused verbatim from build_tissue_markers.py Step 5) ──

def select_markers(norm_matrix: np.ndarray, available: list[str],
                    windows_filt: pd.DataFrame, high: float, low: float,
                    top_n: int, out_dir: Path,
                    blood_veto_mask: np.ndarray,
                    blood_tissues: set,
                    target_fdr: float = 0.05,
                    n_permutations: int = 100,
                    rng_seed: int = 0,
                    min_markers_floor: int = 0,
                    top_n_overrides: dict | None = None,
                    high_pct: float = 0.0,
                    low_pct: float = 0.0) -> tuple[pd.DataFrame, dict]:
    top_n_overrides = top_n_overrides or {}
    tissue_idx = {t: i for i, t in enumerate(available)}
    summary_rows = []
    snr_curves: dict[str, np.ndarray] = {}
    all_marker_idx: set[int] = set()   # union of selected rows = the A/universe row set

    snr_floor = snr_denom_floor(norm_matrix)
    print(f"  SNR denominator floor (10th pct of nonzero signal): {snr_floor:.4g} "
          f"(replaces a fixed 1e-6 epsilon -- see snr_denom_floor docstring)")
    if n_permutations > 0:
        print(f"  Fitting permutation null ({n_permutations} column-shuffle reps, "
              f"target FDR {target_fdr:g}) for the variable top-N cutoff …")
        rng = np.random.default_rng(rng_seed)
        null_snr = permutation_null_snr(
            norm_matrix, available, high, low, blood_veto_mask, blood_tissues,
            n_permutations, rng, snr_floor,
            high_pct=high_pct, low_pct=low_pct)
    else:
        print("  Permutation-null FDR selection: DISABLED (--n-permutations 0) -- "
              "falling back to --high/--low/veto + fixed --top-n.")
        null_snr = {t: np.array([]) for t in available}

    # Row-wise "max of all OTHER columns" computed ONCE for every tissue via
    # the top-2-per-row trick (same one permutation_null_snr already uses
    # above), instead of a fresh np.delete(norm_matrix, j, axis=1).max(axis=1)
    # per tissue -- avoids ~len(available) redundant full-matrix-minus-one-
    # column copies of a multi-million-row array.
    n_cols = norm_matrix.shape[1]
    if n_cols >= 2:
        _argmax_col = np.argmax(norm_matrix, axis=1)
        _top2 = np.partition(norm_matrix, -2, axis=1)[:, -2:]
        _top1, _second = _top2[:, 1], _top2[:, 0]

    # ── Scale-free (per-column quantile) selection ──────────────────────────
    # Step 6 scales every column by an RLM slope fitted at housekeeping anchors.
    # Bedgraph-derived (Roadmap) columns land within ~one order of magnitude of
    # each other (bulk build: 0.60-6.0), but FRAGMENT-derived columns do not:
    # brain_pb came out at 0.071 and the single-cell CNS types at 0.088-2.5
    # against neutrophil's 20.2. A single ABSOLUTE --high/--low pair therefore
    # cannot fit both ingestion paths at once -- on the sparse fragment columns
    # (brain_pb: median 0.000, p99 369) only isolated normalisation spikes clear
    # --high, so the panel fills with non-specific single-bin windows (63.7% of
    # brain_pb's markers were one 200 bp bin; only 4.7% of the genuine bulk
    # brain panel was recovered), and those spikes inflate other_max for every
    # other tissue, collapsing their SNRs (adipose 13.7->6.4, neutrophil
    # 50.3->9.2).
    #
    # With --high-pct/--low-pct each tissue is thresholded at a quantile of its
    # OWN column, which is invariant to that column's scale factor. The
    # specificity test becomes "every OTHER tissue below its own low quantile"
    # rather than "the row max below one global constant" -- a per-column test
    # too, so it cannot be dominated by whichever column happens to be scaled
    # largest.
    #
    # OPT-IN: both must be > 0. Left at 0 the absolute --high/--low path below
    # runs unchanged, so existing bulk builds reproduce byte-for-byte.
    use_pct = high_pct > 0.0 and low_pct > 0.0
    if use_pct:
        # Per column rather than a single np.quantile(axis=0) call: the latter
        # copies the whole 4.3M x n_cols matrix to sort it.
        hi_thr = np.array([np.quantile(norm_matrix[:, k], high_pct)
                           for k in range(n_cols)])
        lo_thr = np.array([np.quantile(norm_matrix[:, k], low_pct)
                           for k in range(n_cols)])
        below_lo = norm_matrix < lo_thr[np.newaxis, :]   # bool, ~1 byte/cell
        n_above  = (~below_lo).sum(axis=1)
        print(f"  Per-column quantile thresholds "
              f"(--high-pct {high_pct:g} / --low-pct {low_pct:g}):")
        for t in available:
            k = tissue_idx[t]
            print(f"    {t:25s}  high={hi_thr[k]:10.3f}  low={lo_thr[k]:10.3f}")

    for tissue in available:
        j = tissue_idx[tissue]
        focal = norm_matrix[:, j]
        if n_cols >= 2:
            other_max = np.where(_argmax_col == j, _second, _top1)
        else:
            other_max = np.zeros_like(focal)

        if use_pct:
            # "every OTHER tissue below its own low quantile": either no column
            # is above its threshold, or the single column above it is the
            # focal one.
            others_ok = (n_above == 0) | ((n_above == 1) & ~below_lo[:, j])
            is_marker = (focal > hi_thr[j]) & others_ok
        else:
            is_marker = (focal > high) & (other_max < low)
        # Blood-open veto: drop candidate windows that sit in blood-accessible
        # chromatin, but only for non-hematopoietic tissues (a blood-open window
        # is a valid marker for the blood lineages themselves).
        n_blood_vetoed = 0
        if tissue not in blood_tissues:
            n_blood_vetoed = int((is_marker & blood_veto_mask).sum())
            is_marker = is_marker & ~blood_veto_mask
        pass_idx = np.where(is_marker)[0]

        snr = focal[pass_idx] / np.maximum(other_max[pass_idx], snr_floor)
        order = np.argsort(snr)[::-1]
        pass_idx, snr = pass_idx[order], snr[order]   # both now SNR-descending
        snr_curves[tissue] = snr.copy()

        # Variable top-N: largest N with permutation-estimated FDR <= target,
        # replacing the old fixed --top-n trim (still available as an optional
        # ceiling below -- e.g. to cap a very large, low-specificity panel).
        n_star, fdr_at_cutoff = fdr_controlled_n(snr, null_snr[tissue],
                                                  n_permutations, target_fdr)
        pass_idx, snr = pass_idx[:n_star], snr[:n_star]

        # Per-tissue --top-n override (see --top-n-overrides): some tissues
        # genuinely support -- and benefit from -- a much larger panel than
        # the rest (e.g. neutrophil's specificity_n_sweep.py curve peaks
        # around a few thousand markers instead of declining from the start
        # like most tissues, because it's the deepest/brightest reference
        # dataset in the panel with a much shallower quality decay). This
        # keeps the DEFAULT --top-n uniform across tissues while allowing a
        # deliberate, documented exception rather than silently forcing
        # every tissue to the same cap.
        effective_top_n = top_n_overrides.get(tissue, top_n)
        if effective_top_n > 0 and len(pass_idx) > effective_top_n:
            pass_idx, snr = pass_idx[:effective_top_n], snr[:effective_top_n]

        # Pool-size exclusion floor (see marker_floor_sweep.py): a tissue whose
        # FDR-controlled panel still comes out below this floor gets EXCLUDED
        # entirely (both its marker rows AND, downstream, its A column via
        # nnls_core.build_reference's zero-owned-rows drop) rather than shipped
        # with a thin, high-variance panel. Reported as `floor_excluded` below
        # rather than silently zeroed, since it changes what tissues exist in
        # the deconvolution, not just how many markers one tissue has.
        floor_excluded = 0 < len(pass_idx) < min_markers_floor
        if floor_excluded:
            pass_idx, snr = pass_idx[:0], snr[:0]

        all_marker_idx.update(pass_idx.tolist())
        marker_df = windows_filt.iloc[pass_idx].copy()
        marker_df["focal_score"] = focal[pass_idx]
        marker_df["snr"] = snr

        n_total = len(marker_df)
        n_tss = (marker_df["type"] == "TSS").sum()
        n_extra = (marker_df["type"] == "EXTRA_TSS").sum()
        n_enh = (marker_df["type"] == "Enhancer").sum()

        bed_path = out_dir / f"{tissue}_markers.bed.gz"
        with gzip.open(bed_path, "wt") as fh:
            if n_total > 0:
                out = marker_df[["chr", "start", "end", "type"]].copy()
                out["score"] = marker_df["focal_score"].round(4)
                out["snr"] = marker_df["snr"].round(4)
                fh.write(out.to_csv(sep="\t", header=False, index=False))

        n_pass_before_trim = len(snr_curves[tissue])
        summary_rows.append({
            "tissue": tissue, "n_pass_thresh": n_pass_before_trim,
            "n_blood_vetoed": n_blood_vetoed,
            "n_star_fdr": n_star, "fdr_at_cutoff": fdr_at_cutoff,
            "n_markers": n_total, "n_TSS": n_tss, "n_EXTRA_TSS": n_extra,
            "n_Enhancer": n_enh,
            "snr_median": float(np.median(snr)) if n_total else 0.0,
            "floor_excluded": floor_excluded,
            "bed": str(bed_path),
        })
        if floor_excluded:
            trim_note = (f"  EXCLUDED: FDR-controlled panel ({n_star} markers) "
                        f"< --min-markers-floor {min_markers_floor}")
        else:
            trim_note = f"  (FDR-cut from {n_pass_before_trim:,} @ FDR={fdr_at_cutoff:.3g})"
            if effective_top_n > 0 and n_star > effective_top_n:
                override_note = (f" [override]" if tissue in top_n_overrides else "")
                trim_note += f"  (further capped to --top-n {effective_top_n}{override_note})"
        print(f"  {tissue:25s}  {n_total:6,} markers{trim_note}  "
              f"(TSS {n_tss:,}  EXTRA_TSS {n_extra:,}  Enhancer {n_enh:,})")

    summary_df = pd.DataFrame(
        summary_rows,
        columns=["tissue", "n_pass_thresh", "n_blood_vetoed", "n_star_fdr",
                 "fdr_at_cutoff", "n_markers", "n_TSS", "n_EXTRA_TSS",
                 "n_Enhancer", "snr_median", "floor_excluded", "bed"])
    marker_idx = np.array(sorted(all_marker_idx), dtype=int)
    return summary_df, snr_curves, marker_idx


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Build Sadeh-background-subtracted H3K27ac tissue markers.")
    ap.add_argument("--cpus", type=int, default=8,
                    help="Parallel replicate workers.")
    ap.add_argument("--high", type=float, default=0.0,
                    help="Focal-tissue threshold (units are now "
                         "background-subtracted counts, not raw pileup — "
                         "re-tune after inspecting the score distribution).")
    ap.add_argument("--low", type=float, default=float("inf"),
                    help="Other-tissue ceiling (see --high note).")
    ap.add_argument("--high-pct", type=float, default=0.0,
                    help="Per-column focal quantile, e.g. 0.99. When >0 (with "
                         "--low-pct) this REPLACES --high: each tissue is "
                         "thresholded at that quantile of its OWN column, so "
                         "selection is invariant to the Step-6 scale factor. "
                         "Required whenever the matrix mixes bedgraph- and "
                         "fragment-derived columns, whose scale factors differ "
                         "by ~300x (brain_pb 0.071 vs neutrophil 20.2). "
                         "0 = keep the absolute --high/--low behaviour.")
    ap.add_argument("--low-pct", type=float, default=0.0,
                    help="Per-column other-tissue quantile, e.g. 0.75. Pairs "
                         "with --high-pct and REPLACES --low: a window "
                         "qualifies only if EVERY other tissue is below its "
                         "own quantile, instead of the row max being below one "
                         "global constant.")
    ap.add_argument("--top-n", type=int, default=0,
                    help="Optional CEILING on markers per tissue (applied after "
                         "the permutation-null FDR cutoff below). 0 = no ceiling. "
                         "Applies uniformly unless overridden per tissue by "
                         "--top-n-overrides.")
    ap.add_argument("--top-n-overrides", type=str, default="",
                    help="Comma-separated tissue=N pairs overriding --top-n for "
                         "specific tissues, e.g. 'neutrophil=300'. Use this for a "
                         "tissue whose own specificity_n_sweep.py curve genuinely "
                         "supports (or needs) a different panel size than the "
                         "rest -- e.g. neutrophil's specificity peaks around a "
                         "few thousand markers instead of declining from the "
                         "start like most tissues, since it's the deepest/"
                         "brightest reference dataset in the panel. Keep this a "
                         "deliberate, documented exception, not a way to re-tune "
                         "every tissue individually.")
    ap.add_argument("--target-fdr", type=float, default=0.05,
                    help="Target FDR for the variable top-N cutoff (only used "
                         "if --n-permutations > 0 -- see that flag's help).")
    ap.add_argument("--n-permutations", type=int, default=0,
                    help="Column-shuffle permutations for the FDR null. "
                         "DEFAULT 0 = DISABLED: this whole-genome marginal-"
                         "preserving permutation null was found to be "
                         "fundamentally unable to distinguish true markers "
                         "from noise on sparse, background-subtracted data -- "
                         "confirmed on synthetic ground truth (a window being "
                         "'zero in every other tissue' isn't rare/distinctive "
                         "when most of the genome is inactive everywhere, so "
                         "relocated real values keep landing on coincidentally "
                         "all-zero windows and scoring just as high as true "
                         "markers, regardless of shuffle scheme or epsilon "
                         "choice -- tried both, neither fixed it). Needs a "
                         "properly redesigned null (e.g. a parametric NB "
                         "significance test reusing the existing Sadeh "
                         "background/dispersion model, BH-corrected, instead "
                         "of empirical resampling) before re-enabling. With "
                         "0, selection falls back to the pre-existing, proven "
                         "--high/--low/veto + fixed --top-n behavior.")
    ap.add_argument("--min-markers-floor", type=int, default=0,
                    help="Exclude a tissue ENTIRELY (drop its column downstream "
                         "too, via nnls_core.build_reference) if its FDR-"
                         "controlled panel comes out below this size. Derive "
                         "this value from marker_floor_sweep.py's variance-"
                         "stabilization sweep, not a guess -- see NNLS_"
                         "deconvolution/scripts/marker_floor_sweep.py. "
                         "0 = no floor (default).")
    ap.add_argument("--window-types", type=str,
                    default="TSS,EXTRA_TSS,Enhancer",
                    help="Comma-separated Windows.csv annotation types to "
                         "treat as signal windows.")
    ap.add_argument("--housekeeping-bed", type=Path,
                    default=HOUSEKEEPING_BED_DEFAULT,
                    help="BED of housekeeping loci used as the "
                         "quantile-normalisation anchor set.")
    ap.add_argument("--bg-hierarchy", choices=["regional", "chromosome", "global"],
                    default="regional",
                    help="Sadeh background-model depth.")
    ap.add_argument("--blood-peaks-broad", type=Path, nargs="+",
                    default=[PBMC_BROAD_BED_DEFAULT, MEGAKARYOCYTE_BROAD_BED_DEFAULT],
                    help="Blood-open peak BEDs used AS-IS (broad calls, e.g. PBMC).")
    ap.add_argument("--blood-peaks-narrow", type=Path, nargs="+",
                    default=[NEUTROPHIL_PEAKS_DEFAULT],
                    help="Blood-open narrow-peak BEDs, slopped by --blood-slop-bp "
                         "to approximate broad coverage (e.g. neutrophil).")
    ap.add_argument("--blood-slop-bp", type=int, default=200,
                    help="Bp to extend each narrow blood peak on both sides.")
    ap.add_argument("--blood-veto-frac", type=float, default=0.5,
                    help="Drop a non-blood candidate window when >= this fraction "
                         "of it is covered by the blood-open mask.")
    ap.add_argument("--chrom-sizes", type=Path, default=CHROM_SIZES_DEFAULT,
                    help="Genome file for bedtools slop (hg19).")
    ap.add_argument("--no-blood-veto", action="store_true",
                    help="Disable the blood-open veto (reproduce old behaviour).")
    ap.add_argument("--fasta", type=Path, default=FASTA_DEFAULT,
                    help="hg19 genome FASTA for GC-content computation "
                         "(bedtools nuc).")
    ap.add_argument("--no-gc-correction", action="store_true",
                    help="Disable the GC-bias residual correction (reproduce "
                         "pre-GC-correction behaviour).")
    ap.add_argument("--cache-dir", type=Path, default=None,
                    help="Fragment/count cache directory "
                         "(default: <out>/_cache).")
    ap.add_argument("--out", type=Path,
                    default=CFCHROM_DIR / "tissue_markers_sadeh",
                    help="Output directory.")
    args = ap.parse_args()

    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.cache_dir or (out_dir / "_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)

    window_types = [t.strip() for t in args.window_types.split(",") if t.strip()]
    bad_types = set(window_types) - SELECTABLE_TYPES
    if bad_types:
        raise ValueError(f"--window-types must be a subset of "
                          f"{sorted(SELECTABLE_TYPES)}, got {sorted(bad_types)}")

    # ── Step 1: windows ───────────────────────────────────────────────────
    print(f"\n[Step 1] Loading signal windows ({window_types}) from Windows.csv …")
    sig_df = load_signal_windows(window_types)
    print(f"  {len(sig_df):,} unique signal windows "
          f"(TSS {(sig_df['type']=='TSS').sum():,}, "
          f"EXTRA_TSS {(sig_df['type']=='EXTRA_TSS').sum():,}, "
          f"Enhancer {(sig_df['type']=='Enhancer').sum():,})")

    print("\n[Step 1b] Loading autosomal BACKGROUND tiles …")
    bg_df = load_background_tiles()
    print(f"  {len(bg_df):,} background tiles")

    # Signal and background windows are disjoint by construction — Windows.csv
    # partitions the genome into non-overlapping TSS/EXTRA_TSS/Enhancer/
    # FLANKING/BACKGROUND tiles, so no subtract step is needed here (unlike
    # module4's curated off-target BED, which does overlap signature regions).

    # ── Step 2: replicate resolution ──────────────────────────────────────
    print("\n[Step 2] Resolving tissue replicate files …")
    tissue_replicates = build_tissue_replicates()
    for tissue, reps in sorted(tissue_replicates.items()):
        print(f"  {tissue:20s}  {len(reps)} replicate(s): "
              f"{[r['id'] for r in reps]}")
    if not tissue_replicates:
        print("ERROR: no tissue replicate files found.")
        return

    # ── Step 3: combined counting BED (signal + background, shared) ───────
    wt_sig = "-".join(sorted(window_types))
    combined_bed = cache_dir / f"combined_windows.{wt_sig}.bed"
    if not combined_bed.exists():
        print("\n[Step 3] Building combined signal+background counting BED …")
        tmp_unsorted = tempfile.NamedTemporaryFile(
            suffix=".bed", delete=False, mode="w")
        for _, r in sig_df.iterrows():
            tmp_unsorted.write(
                f"{r['chr']}\t{int(r['start'])}\t{int(r['end'])}\t{r['region_id']}\n")
        for _, r in bg_df.iterrows():
            tmp_unsorted.write(
                f"{r['chr']}\t{int(r['start'])}\t{int(r['end'])}\t{r['region_id']}\n")
        tmp_unsorted.close()
        run(f"sort -k1,1 -k2,2n {tmp_unsorted.name} > {combined_bed}")
        Path(tmp_unsorted.name).unlink()
    print(f"  Combined BED → {combined_bed}")

    # GC content is a property of the genome, not the replicate — computed
    # ONCE for every window (signal + background) and reused by every
    # replicate's own per-replicate GC-correction fit (see process_replicate).
    gc_correction = not args.no_gc_correction
    if gc_correction:
        print(f"\n[Step 3b] Computing GC content ({args.fasta.name}) …")
        gc_map = compute_gc_content(str(combined_bed), str(args.fasta))
        sig_df = sig_df.copy()
        sig_df["gc"] = sig_df["region_id"].map(gc_map)
        bg_gc = bg_df["region_id"].map(gc_map).to_numpy()
        print(f"  GC computed for {len(gc_map):,} windows")
    else:
        print("\n[Step 3b] GC-bias correction disabled (--no-gc-correction).")
        bg_gc = np.zeros(len(bg_df))

    bg_arrays = {
        "region_id": bg_df["region_id"].to_numpy(),
        "length_kb": bg_df["length_kb"].to_numpy(),
        "chr": bg_df["chr"].to_numpy(),
        "region_2mb": bg_df["region_2mb"].to_numpy(),
        "gc": bg_gc,
    }

    # ── Step 4: per-replicate counting + background subtraction (parallel) ─
    print(f"\n[Step 4] Counting + background-subtracting "
          f"({sum(len(r) for r in tissue_replicates.values())} replicates, "
          f"{args.cpus} workers) …")
    flat_reps = [
        (rep, tissue)
        for tissue, reps in tissue_replicates.items()
        for rep in reps
    ]
    per_tissue_signals: dict[str, list[np.ndarray]] = {t: [] for t in tissue_replicates}
    ref_raw_counts: list[np.ndarray] = []   # per-replicate raw signal-window counts
    with ThreadPoolExecutor(max_workers=args.cpus) as pool:
        futures = {
            pool.submit(process_replicate, rep, tissue, combined_bed, wt_sig,
                        sig_df, bg_arrays, args.bg_hierarchy, cache_dir,
                        gc_correction): (rep, tissue)
            for rep, tissue in flat_reps
        }
        for fut in as_completed(futures):
            tissue, _rep_id, signal, raw_counts = fut.result()
            per_tissue_signals[tissue].append(signal)
            ref_raw_counts.append(raw_counts)

    # ── Step 4b: mean-variance diagnostic on raw reference window counts ────
    #   Note: variance is taken across ALL replicates of ALL tissues, so it
    #   mixes technical noise with genuine cross-tissue biological signal — read
    #   the fitted alpha as an upper bound on pure count overdispersion.
    if len(ref_raw_counts) >= 2:
        mean_variance_plot(
            np.column_stack(ref_raw_counts),
            out_dir / "mean_variance_reference.png",
            f"Window-count mean–variance across "
            f"{len(ref_raw_counts)} reference replicates",
            "signal windows, across replicates")

    # ── Step 4c: GC-correction diagnostics (only when correction is active) ─
    # Proves the GC correction applied in process_replicate is warranted (does
    # GC bias exist?) and non-harmful (does correcting it move a GC-matched
    # fake-signature null toward 0?). Reference-side analogue of module4's
    # within-sample GC checks. Skipped under --no-gc-correction (no gc column).
    if gc_correction:
        print("\n[Step 4c] GC-bias correction diagnostics "
              "(residual-vs-GC + fake-signature null) …")
        try:
            gc_correction_diagnostics(
                flat_reps, combined_bed, wt_sig, sig_df, bg_arrays,
                args.bg_hierarchy, cache_dir, out_dir)
        except Exception as _gc_diag_exc:
            print(f"  [gc-diag] skipped due to error: {_gc_diag_exc}")

    # ── Step 5: average replicates -> raw (background-subtracted) matrix ───
    print("\n[Step 5] Averaging replicates per tissue …")
    available = sorted(t for t, sigs in per_tissue_signals.items() if sigs)
    raw_matrix = np.column_stack([
        np.mean(np.stack(per_tissue_signals[t]), axis=0) for t in available
    ]).astype(np.float64)
    print(f"  {len(available)} tissues available: {available}")

    # ── Step 6: housekeeping-anchored quantile normalisation ───────────────
    print(f"\n[Step 6] Quantile-normalising against {args.housekeeping_bed} …")
    tmp_win_bed = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    for idx, r in sig_df.iterrows():
        tmp_win_bed.write(f"{r['chr']}\t{r['start']}\t{r['end']}\t{idx}\n")
    tmp_win_bed.close()
    tmp_win_sorted = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    run(f"sort -k1,1 -k2,2n {tmp_win_bed.name} > {tmp_win_sorted}")
    Path(tmp_win_bed.name).unlink()
    hk_out = run(
        f"bedtools intersect -a {tmp_win_sorted} -b {args.housekeeping_bed} -u"
    )
    Path(tmp_win_sorted).unlink()
    anchor_idx = [int(line.split("\t")[3]) for line in hk_out.splitlines()]
    anchor_mask = np.zeros(len(sig_df), dtype=bool)
    anchor_mask[anchor_idx] = True
    n_anchor = int(anchor_mask.sum())
    print(f"  {n_anchor:,} signal windows overlap housekeeping loci")

    if n_anchor < 10:
        print("  WARNING: fewer than 10 anchor windows — skipping "
              "normalisation, using raw background-subtracted matrix.")
        norm_matrix = raw_matrix.copy()
        scale_factors = np.ones(raw_matrix.shape[1])
    else:
        sub_raw = raw_matrix[anchor_mask]
        sub_qq = quantile_normalize(sub_raw)
        scale_factors = np.array([
            fit_rlm_slope(sub_raw[:, j], sub_qq[:, j])
            for j in range(raw_matrix.shape[1])
        ])
        norm_matrix = raw_matrix * scale_factors[np.newaxis, :]

    print("  Scale factors:")
    for tissue, sf in zip(available, scale_factors):
        print(f"    {tissue:25s}  {sf:.6f}")

    norm_df = pd.DataFrame(norm_matrix, columns=available)
    norm_df.insert(0, "type", sig_df["type"].values)
    norm_df.insert(0, "end", sig_df["end"].values)
    norm_df.insert(0, "start", sig_df["start"].values)
    norm_df.insert(0, "chr", sig_df["chr"].values)
    norm_tsv = out_dir / "normalised_counts.tsv.gz"
    norm_df.to_csv(norm_tsv, sep="\t", index=False, compression="gzip")
    print(f"  Normalised counts → {norm_tsv}")

    # ── Step 6b: blood-open veto mask ─────────────────────────────────────
    if args.no_blood_veto:
        blood_veto_mask = np.zeros(len(sig_df), dtype=bool)
        print("\n[Step 6b] Blood-open veto disabled (--no-blood-veto).")
    else:
        print(f"\n[Step 6b] Building blood-open veto mask "
              f"(PBMC broad + neutrophil narrow slop±{args.blood_slop_bp} bp, "
              f"veto ≥ {args.blood_veto_frac} window overlap) …")
        blood_veto_mask = build_blood_veto_mask(
            sig_df, args.blood_peaks_broad, args.blood_peaks_narrow,
            args.blood_slop_bp, args.chrom_sizes, args.blood_veto_frac,
            cache_dir)
        vetoed_tissues = sorted(set(available) - BLOOD_TISSUES)
        print(f"  {int(blood_veto_mask.sum()):,} / {len(sig_df):,} signal "
              f"windows flagged blood-open")
        print(f"  vetoed for non-blood tissues: {vetoed_tissues}")
        print(f"  exempt (hematopoietic): "
              f"{sorted(set(available) & BLOOD_TISSUES)}")

    # ── Step 7: marker selection ─────────────────────────────────────────
    top_n_overrides = {}
    for pair in args.top_n_overrides.split(","):
        pair = pair.strip()
        if not pair:
            continue
        tissue, _, n_str = pair.partition("=")
        tissue = tissue.strip()
        if not n_str or tissue not in available:
            raise ValueError(
                f"--top-n-overrides entry {pair!r} is malformed or names a "
                f"tissue not in {available}")
        top_n_overrides[tissue] = int(n_str)

    fdr_desc = (f"target FDR {args.target_fdr:g} via {args.n_permutations} permutations"
                if args.n_permutations > 0 else "FDR selection DISABLED")
    _use_pct = args.high_pct > 0.0 and args.low_pct > 0.0
    _thr_desc = (f"focal > per-column q{args.high_pct:g}, "
                 f"all others < their own q{args.low_pct:g}"
                 if _use_pct else
                 f"focal > {args.high}, all others < {args.low}")
    print(f"\n[Step 7] Marker selection "
          f"({_thr_desc}, {fdr_desc}"
          + (f", capped at top-{args.top_n}" if args.top_n > 0 else "")
          + (f", overrides {top_n_overrides}" if top_n_overrides else "")
          + (f", floor {args.min_markers_floor}" if args.min_markers_floor > 0 else "")
          + ") …")
    summary_df, snr_curves, marker_idx = select_markers(
        norm_matrix, available, sig_df, args.high, args.low, args.top_n, out_dir,
        blood_veto_mask, BLOOD_TISSUES,
        target_fdr=args.target_fdr, n_permutations=args.n_permutations,
        min_markers_floor=args.min_markers_floor,
        top_n_overrides=top_n_overrides,
        high_pct=args.high_pct, low_pct=args.low_pct)

    n_floor_excluded = int(summary_df["floor_excluded"].sum())
    if n_floor_excluded:
        excluded_tissues = summary_df.loc[summary_df["floor_excluded"], "tissue"].tolist()
        print(f"  {n_floor_excluded} tissue(s) EXCLUDED by --min-markers-floor "
              f"{args.min_markers_floor}: {excluded_tissues}")

    qc_dir = out_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    qc_path = qc_dir / "permutation_null_summary.tsv"
    summary_df.to_csv(qc_path, sep="\t", index=False)
    print(f"  Permutation-null FDR summary → {qc_path}")

    # ── Step 7b: tissue-tissue correlation heatmap on the marker rows ─────
    # Computed on the union of selected markers (= the universe / A row set) so
    # it reflects the columns NNLS will actually try to separate.
    print("\n[Step 7b] Tissue correlation heatmap on the marker set …")
    if marker_idx.size >= 2:
        correlation_heatmap(
            norm_matrix[marker_idx, :], available,
            out_dir / "marker_correlation_heatmap.png",
            f"Tissue signature correlation on {marker_idx.size:,} marker windows")
    else:
        print("  too few markers selected; skipping correlation heatmap.")

    # ── Step 8: SNR rank-curve plot ───────────────────────────────────────
    plot_path = out_dir / "snr_rank_curves.png"
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.get_cmap("tab20", max(len(snr_curves), 1))
    for idx, (tissue, curve) in enumerate(sorted(snr_curves.items())):
        ranks = np.arange(1, len(curve) + 1)
        ax.plot(ranks, curve, lw=1.2, label=tissue, color=cmap(idx))
    if args.top_n > 0:
        ax.axvline(args.top_n, color="black", ls="--", lw=1.2,
                   label=f"top-N = {args.top_n}")
    ax.set_xlabel("Rank (by SNR, descending)")
    ax.set_xlim(xmin=0, xmax=2000)
    ax.set_ylabel("SNR  [focal / max(others)]")
    ax.set_title("SNR rank curves per tissue — Sadeh background-subtracted")
    ax.set_yscale("log")
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  SNR rank curves → {plot_path}")

    # ── Step 9: summary ────────────────────────────────────────────────────
    summary_tsv = out_dir / "marker_summary.tsv"
    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    print(f"\n[Summary] → {summary_tsv}")
    print(summary_df.to_string(index=False))
    print(f"\nDone. Output directory: {out_dir}")


if __name__ == "__main__":
    main()
