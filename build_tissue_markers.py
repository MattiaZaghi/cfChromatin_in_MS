#!/usr/bin/env python3
"""
build_tissue_markers.py  —  H3K27ac tissue-specific marker BED files.

Pipeline
--------
1.  Read Windows.csv; keep rows with annotation type ∈ {TSS, EXTRA_TSS, Enhancer}.
    Collapse duplicate coordinates (same chr:start-end) to the highest-priority
    type: TSS > EXTRA_TSS > Enhancer.
2.  For multi-sample tissues (t_cell), build a mean-signal bedgraph from
    individual sample bedgraphs (bedtools unionbedg + per-position mean).
3.  For each tissue, sum bedgraph signal over every window (bedtools map -o sum).
4.  QQ-normalise the resulting windows × tissues matrix:
      (a) Constitutive windows = those where every tissue has signal > QQ_FLOOR.
      (b) Quantile-normalise across tissues for those windows.
      (c) Fit a per-tissue linear scale factor via robust regression
          (Huber M-estimator; OLS fallback) on QQ-normalised ~ raw.
      (d) Apply scale factors to the full matrix.
5.  For each tissue: windows with focal score > HIGH_THR and every other tissue
    score < LOW_THR are tissue-specific markers.
6.  Write one gzipped BED per tissue + a summary TSV + console report.

Notes
-----
- Bedgraph normalisation: signal = coverage_per_bp × 10000 / total_reads
  (RC = Read-Count scaled, matching create_brain_bedgraph.py).
- Missing bedgraph files → tissue is skipped with a warning.
- Breast, Prostate, Skeletal-muscle-F, Megakaryocyte are marked TODO;
  add their bedgraph filenames to TISSUES once the files are created.

Usage
-----
    python build_tissue_markers.py [--cpus N] [--high THR] [--low THR] [--top-n N]

    --cpus N    Parallel bedtools map workers (default: 8).
    --high THR  Focal-tissue signal threshold (default: 35).
    --low THR   Other-tissue ceiling (default: 15).
    --floor F   Minimum signal to count as constitutive for QQ (default: 1e-4).
    --top-n N   Keep top-N markers per tissue ranked by focal/max(others) SNR.
                0 = keep all passing markers (default: 0).
    --out DIR   Output directory (default: ./tissue_markers).
"""

import argparse
import gzip
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

# ── Paths ───────────────────────────────────────────────────────────────────
BASE        = Path("/cfs/klemming/projects/supr/uppstore2017150/anita")
WINDOWS_CSV = BASE / "Sadeh_markers/Windows.csv"
BGDIR       = BASE / "epiroadmap/bedgraph"

# ── Annotation type priority (higher = preferred when deduplicating) ─────────
ANN_PRIORITY = {"Enhancer": 0, "EXTRA_TSS": 1, "TSS": 2}
KEEP_TYPES   = set(ANN_PRIORITY)

# ── Tissue → bedgraph file(s) mapping ────────────────────────────────────────
# Lists with one file → used directly.
# Lists with multiple files → merged on-the-fly (mean signal).
# Comment out or remove entries whose bedgraphs don't exist yet.
TISSUES = {
    "adipose":          ["Adipose_E063_RC.bedgraph"],
    "b_cell":           ["B_naive_E032_RC.bedgraph"],
    "brain":            ["Brain_RC.bedgraph"],
    "breast":           ["Breast_RC.bedgraph"],
    "epithelial":       ["Epithelial_E119_RC.bedgraph"],
    "digestive":        ["digestive_merged_RC.bedgraph"],  # pre-merged Colon+Rectum+Stomach
    "heart":            ["heart_merged_RC.bedgraph"],      # pre-merged Atrium+Ventricle+Aorta
    "t_cell":           [
        "CD4_Th1_E041_RC.bedgraph",
        "CD4_Th17_E042_RC.bedgraph",
        "CD4_Treg_E044_RC.bedgraph",
        "CD4_T_naive_E038_RC.bedgraph",
        "CD4_T_memory_E037_RC.bedgraph",
        "CD8_T_naive_E047_RC.bedgraph",
        "CD8_T_memory_E048_RC.bedgraph",
    ],
    "liver":            ["Liver_E066_RC.bedgraph"],
    "lung":             ["Lung_E128_RC.bedgraph"],
    "pancreas":         ["pancreas_merged_RC.bedgraph"],   # pre-merged Pancreas+Islet
    "prostate":         ["Prostate_RC.bedgraph"],
    "ovary":            ["Ovary_E097_RC.bedgraph"],
    "nk":               ["NK_E046_RC.bedgraph"],
    "neutrophil":       ["Neutrophil_ENCFF067WAA.bedGraph"],
    "monocyte":         ["Monocyte_CD14_E029_RC.bedgraph"],
    "skin":             ["Skin_E127_RC.bedgraph"],
    "placenta":         ["Placenta_E091_RC.bedgraph"],
    "gi_smooth_muscle": ["GI_Smooth_Muscle_E111_RC.bedgraph"],
    "skeletal_muscle_f":["Sk_muscle_F_E108_RC.bedgraph"]
    # "megakaryocyte":    ["Megakaryocyte_RC.bedgraph"],    # BLUEPRINT TODO
}

# narrowpeak files
REFERENCE_DATA = {
    # ── Roadmap Epigenomics — immune cells (hg19 native) ─────────────────
    # Functional tags ("immune", "bcell") must appear first — they drive the
    # extract_exclusive_subset grep logic in Step 4.
    # Descriptive tags from Supplemental Table 12 Groups column follow.
    # EID listed as comment after the URL.
    "B_naive":         {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E032-H3K27ac.narrowPeak.gz",   # E032 BLD.CD19.PPC
                        "tags": ["bcell", "immune", "B-Cells", "Lymphocytes", "Leukocytes"]},
    "CD4_Th1":         {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E041-H3K27ac.narrowPeak.gz",   # E041 BLD.CD4.CD25M.IL17M.PL.TPC
                        "tags": ["immune", "T-Helper-Cells", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "CD4_Th17":        {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E042-H3K27ac.narrowPeak.gz",   # E042 BLD.CD4.CD25M.IL17P.PL.TPC
                        "tags": ["immune", "T-Helper-Cells", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "CD4_Treg":        {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E044-H3K27ac.narrowPeak.gz",   # E044 BLD.CD4.CD25.CD127M.TREGPC
                        "tags": ["immune", "Treg-Cells", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "CD4_T_naive":     {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E038-H3K27ac.narrowPeak.gz",   # E038 BLD.CD4.NPC
                        "tags": ["immune", "T-Helper-Cells", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "CD4_T_memory":    {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E037-H3K27ac.narrowPeak.gz",   # E037 BLD.CD4.MPC
                        "tags": ["immune", "T-Helper-Cells", "T-Memory-Cells", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "CD8_T_naive":     {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E047-H3K27ac.narrowPeak.gz",   # E047 BLD.CD8.NPC
                        "tags": ["immune", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "CD8_T_memory":    {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E048-H3K27ac.narrowPeak.gz",   # E048 BLD.CD8.MPC
                        "tags": ["immune", "T-Memory-Cells", "T-Cells", "Lymphocytes", "Leukocytes"]},
    "Monocyte_CD14":   {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E029-H3K27ac.narrowPeak.gz",   # E029 BLD.CD14.PC
                        "tags": ["immune", "Monocytes", "Leukocytes"]},
    "NK":              {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E046-H3K27ac.narrowPeak.gz",   # E046 BLD.CD56.PC
                        "tags": ["immune", "NK", "Lymphocytes", "Leukocytes"]},
    "PBMC":            {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E062-H3K27ac.narrowPeak.gz",   # E062 BLD.PER.MONUC.PC
                        "tags": ["immune", "PBMC", "Leukocytes"]},
    "Neutrophil":      {"url": "https://www.encodeproject.org/files/ENCFF457CTA/@@download/ENCFF457CTA.bed.gz",
                        "tags": ["immune", "Neutrophils", "Leukocytes"]},
    "Brain":            {"url": "/cfs/klemming/projects/supr/uppstore2017150/anita/epiroadmap/merge_narrowpeaks_brain/Brain.bed.gz",   # E062 BLD.PER.MONUC.PC
                        "tags": ["cns"]},
    # ── Roadmap Epigenomics — TODO placeholders (BLUEPRINT replacements) ──
    "B_memory":        {"url": "TODO_BLUEPRINT_B_memory_H3K27ac_peaks.bed.gz",
                        "tags": ["bcell", "immune", "B-Cells", "Lymphocytes", "Leukocytes"]},
    "Monocyte_CD16":   {"url": "TODO_BLUEPRINT_CD16mono_H3K27ac_peaks.bed.gz",
                        "tags": ["immune", "Monocytes", "Leukocytes"]},
    "Megakaryocyte":   {"url": "TODO_BLUEPRINT_Megakaryocyte_H3K27ac_peaks.bed.gz",
                        "tags": ["immune", "Megakaryocytes"]},
    # ── Roadmap Epigenomics — non-immune tissues (hg19 native) ───────────
    "Adipose":         {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E063-H3K27ac.narrowPeak.gz",   # E063 FAT.ADIP.NUC
                        "tags": ["Adipose"]},
    "Atrium_Heart":    {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E104-H3K27ac.narrowPeak.gz",   # E104 HRT.ATR.R
                        "tags": ["Atrium", "Heart"]},
    "Ventricle_Heart": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E105-H3K27ac.narrowPeak.gz",   # E105 HRT.VNT.R
                        "tags": ["Ventricle", "Heart"]},
    "Aorta_Heart":     {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E065-H3K27ac.narrowPeak.gz",   # E065 VAS.AOR
                        "tags": ["Aorta", "Heart"]},
    "GI_Smooth_Muscle":{"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E111-H3K27ac.narrowPeak.gz",   # E111 GI.STMC.MUS
                        "tags": ["GI_Sm_Muscle", "GI Sm. Muscle"]},
    "Colon_Digestive": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E106-H3K27ac.narrowPeak.gz",   # E106 GI.CLN.SIG
                        "tags": ["Colon", "Digestive"]},
    "Rectum_GI_Mucosa":{"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E102-H3K27ac.narrowPeak.gz",   # E102 GI.RECT.MUC.31
                        "tags": ["Rectum", "GI Mucosa", "Digestive"]},
    "Stomach_Digestive":{"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E094-H3K27ac.narrowPeak.gz",  # E094 GI.STMC.GAST
                        "tags": ["Stomach", "Digestive"]},
    "Placenta":        {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E091-H3K27ac.narrowPeak.gz",   # E091 PLCNT.FET
                        "tags": ["Placenta"]},
    "Ovary":           {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E097-H3K27ac.narrowPeak.gz",   # E097 OVRY
                        "tags": ["Ovary"]},
    "Pancreas_Islet":  {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E087-H3K27ac.narrowPeak.gz",   # E087 PANC.ISLT
                        "tags": ["Pancreas Islet", "Pancreas"]},
    "Liver":           {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E066-H3K27ac.narrowPeak.gz",   # E066 LIV.ADLT
                        "tags": ["Liver"]},
    "Pancreas":        {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E098-H3K27ac.narrowPeak.gz",   # E098 PANC
                        "tags": ["Pancreas"]},
    "Lung":            {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E128-H3K27ac.narrowPeak.gz",   # E128 LNG.NHLF
                        "tags": ["Lung"]},
    "Epithelial":      {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E119-H3K27ac.narrowPeak.gz",   # E119 BRST.HMEC
                        "tags": ["Epithelial"]},
    "Vasculary":       {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E122-H3K27ac.narrowPeak.gz",   # E122 VAS.HUVEC
                        "tags": ["Vasculary"]},
    "Skin":            {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E127-H3K27ac.narrowPeak.gz",   # E127 SKIN.NHEK
                        "tags": ["Skin"]},
    "Breast_51":          {"url": "https://www.encodeproject.org/files/ENCFF569TYE/@@download/ENCFF569TYE.bed.gz",   
                        "tags": ["Breast"]},
    "Breast_53":          {"url": "https://www.encodeproject.org/files/ENCFF154XFN/@@download/ENCFF154XFN.bed.gz",
                           "tags": ["Breast"]},
    "Prostate_37":         {"url": "https://www.encodeproject.org/files/ENCFF002NFG/@@download/ENCFF002NFG.bed.gz",   
                        "tags": ["Prostate"]},
    "Prostate_54":         {"url": "https://www.encodeproject.org/files/ENCFF063VHN/@@download/ENCFF063VHN.bed.gz",
                           "tags": ["Prostate"]},
    "Sk_muscle_F":      {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E108-H3K27ac.narrowPeak.gz",   # E120 SKIN.FIBR
                        "tags": ["Sk_muscle"]},
    #"Psoas_muscle":      {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E100-H3K27ac.narrowPeak.gz",   # E109 PSOA.MYO
                        #"tags": ["Psoas_muscle"]}
    # ── Anita single-cell (CNS) ───────────────────────────────────────────────
    #"Astrocyte":{"url": "/cfs/klemming/projects/supr/uppstore2017150/anita/all_marker_bins_consistent_hg19/AST_brain_wb.bed.gz",
                       #"tags": ["cns"]},
    #"Microglia":         {"url": "/cfs/klemming/projects/supr/uppstore2017150/anita/all_marker_bins_consistent_hg19/MIGL_brain_wb.bed.gz",
                       #"tags": ["cns"]},
    #"Oligodendrocyte":      {"url": "/cfs/klemming/projects/supr/uppstore2017150/anita/all_marker_bins_consistent_hg19/MOL_brain_wb.bed.gz",
                       #"tags": ["cns"]},
    #"Neuron":      {"url": "/cfs/klemming/projects/supr/uppstore2017150/anita/all_marker_bins_consistent_hg19/NEU_brain_wb.bed.gz",
                       #"tags": ["cns"]},
    #"OPC":      {"url": "/cfs/klemming/projects/supr/uppstore2017150/anita/all_marker_bins_consistent_hg19/OPC_brain_wb.bed.gz",
                       #"tags": ["cns"]}
    #"VLMC":      {"url": "/date/gcb/gcb_MZ/MS-cfEpiMap/all_markers/VLMC.bed.gz",
                       #"tags": ["cns"]}
}



# ── Helpers ──────────────────────────────────────────────────────────────────

def run(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def merge_bedgraphs(bg_paths: list[Path], out: Path) -> bool:
    """Build mean-signal bedgraph from a list of sorted bedgraphs.
    Matches the logic in create_brain_bedgraph.py.
    """
    if out.exists():
        print(f"    exists: {out.name}")
        return True
    missing = [p for p in bg_paths if not p.exists()]
    if missing:
        print(f"    SKIP — missing: {[p.name for p in missing]}")
        return False
    n = len(bg_paths)
    bg_list = " ".join(str(p) for p in bg_paths)
    print(f"    merging {n} bedgraphs → {out.name} …")
    try:
        run(
            f"bedtools unionbedg -i {bg_list} -filler 0 "
            f"| awk 'BEGIN{{OFS=\"\\t\"}}"
            f"{{s=0; for(i=4;i<=NF;i++) s+=$i; print $1,$2,$3,s/{n}}}' "
            f"> {out}"
        )
        return True
    except subprocess.CalledProcessError:
        if out.exists():
            out.unlink()
        return False


def bedtools_map(windows_bed: Path, bedgraph: Path, n_windows: int) -> np.ndarray:
    """Return a float array of summed bedgraph signal per window row."""
    result = subprocess.run(
        f"bedtools map -a {windows_bed} -b {bedgraph} -c 4 -o sum -null 0",
        shell=True, check=True, capture_output=True, text=True,
    )
    lines = result.stdout.splitlines()
    if len(lines) != n_windows:
        raise RuntimeError(
            f"bedtools map returned {len(lines)} lines for {bedgraph.name}, "
            f"expected {n_windows}. stderr: {result.stderr[:300]}"
        )
    return np.fromiter(
        (float(line.rsplit("\t", 1)[-1]) for line in lines),
        dtype=float, count=n_windows,
    )


# ── QQ normalization ─────────────────────────────────────────────────────────

def quantile_normalize(X: np.ndarray) -> np.ndarray:
    """Quantile-normalise columns of X (n_windows × n_tissues).
    Each column gets the same rank-ordered mean values as the pooled distribution.
    """
    sorted_X   = np.sort(X, axis=0)
    row_means  = sorted_X.mean(axis=1)
    result     = np.empty_like(X, dtype=float)
    for j in range(X.shape[1]):
        rank_order = np.argsort(X[:, j], kind="stable")
        result[rank_order, j] = row_means
    return result


def fit_rlm_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Slope of robust regression y ~ x with no intercept (Huber M-estimator)."""
    x_col = x.reshape(-1, 1)
    try:
        model  = sm.RLM(y, x_col, M=sm.robust.norms.HuberT())
        return float(model.fit(disp=False).params[0])
    except Exception:
        # OLS fallback: slope = (x·y) / (x·x)
        denom = float(np.dot(x, x))
        return float(np.dot(x, y)) / denom if denom > 0 else 1.0


def qq_normalize(raw: np.ndarray, floor: float) -> tuple[np.ndarray, np.ndarray]:
    """QQ-normalise a windows × tissues matrix.

    Returns (normalised_matrix, scale_factors).
    """
    # (a) constitutive windows: signal > floor in every tissue
    constitutive = np.all(raw > floor, axis=1)
    n_const = int(constitutive.sum())
    print(f"  constitutive windows (signal > {floor} in all tissues): "
          f"{n_const:,} / {raw.shape[0]:,}")

    if n_const < 200:
        print("  WARNING: too few constitutive windows — skipping QQ, "
              "returning raw matrix with unit scale factors.")
        return raw.copy(), np.ones(raw.shape[1])

    # (b) quantile-normalise constitutive sub-matrix
    sub_raw = raw[constitutive]          # (n_const, n_tissues)
    sub_qq  = quantile_normalize(sub_raw)

    # (c) per-tissue robust scale factor: QQ ~ raw * scale
    n_tissues = raw.shape[1]
    scale_factors = np.array([
        fit_rlm_slope(sub_raw[:, j], sub_qq[:, j])
        for j in range(n_tissues)
    ])

    # (d) apply to full matrix
    normalised = raw * scale_factors[np.newaxis, :]
    return normalised, scale_factors


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Build tissue-specific H3K27ac marker BEDs.")
    ap.add_argument("--cpus",  type=int,   default=8,    help="Parallel bedtools workers.")
    ap.add_argument("--high",  type=float, default=35.0, help="Focal-tissue threshold.")
    ap.add_argument("--low",   type=float, default=15.0, help="Other-tissue ceiling.")
    ap.add_argument("--floor", type=float, default=1e-4, help="QQ constitutive floor.")
    ap.add_argument("--top-n", type=int,   default=0,
                    help="Keep top-N markers per tissue by focal/max(others) SNR. "
                         "0 = keep all passing markers.")
    ap.add_argument("--out",   type=Path,
                    default=Path("/cfs/klemming/projects/supr/uppstore2017150/anita"
                                 "/cfChromatin_in_MS/tissue_markers"),
                    help="Output directory.")
    args = ap.parse_args()
    OUT_DIR = args.out
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load and filter Windows.csv ──────────────────────────────────
    print("\n[Step 1] Loading Windows.csv …")
    # Columns: chr, start, end, gene, tissue_annotation, type
    windows_raw = pd.read_csv(
        WINDOWS_CSV, header=None,
        names=["chr", "start", "end", "gene", "tissue_ann", "type"],
        dtype={"chr": str, "start": int, "end": int,
               "gene": str, "tissue_ann": str, "type": str},
    )
    keep = windows_raw["type"].isin(KEEP_TYPES)
    windows_filt = windows_raw[keep].copy()
    print(f"  rows in KEEP_TYPES: {len(windows_filt):,}  "
          f"(TSS: {(windows_filt['type']=='TSS').sum():,}, "
          f"EXTRA_TSS: {(windows_filt['type']=='EXTRA_TSS').sum():,}, "
          f"Enhancer: {(windows_filt['type']=='Enhancer').sum():,})")

    # Deduplicate by coordinates: keep highest-priority annotation type.
    windows_filt["_pri"] = windows_filt["type"].map(ANN_PRIORITY)
    windows_filt = (
        windows_filt
        .sort_values("_pri", ascending=False)
        .drop_duplicates(subset=["chr", "start", "end"], keep="first")
        .drop(columns="_pri")
        .sort_values(["chr", "start"])
        .reset_index(drop=True)
    )
    print(f"  unique windows after dedup: {len(windows_filt):,}")

    # Write sorted BED for bedtools (vectorised — faster than iterrows)
    tmp_windows = tempfile.NamedTemporaryFile(suffix=".bed", delete=False)
    tmp_windows.close()
    windows_bed = Path(tmp_windows.name)
    windows_filt[["chr", "start", "end", "type"]].to_csv(
        windows_bed, sep="\t", header=False, index=False
    )

    # ── Step 2: Resolve / merge tissue bedgraphs ──────────────────────────────
    print("\n[Step 2] Resolving tissue bedgraphs …")
    tissue_bedgraphs: dict[str, Path] = {}
    for tissue, files in TISSUES.items():
        paths = [BGDIR / f for f in files]
        if len(paths) == 1:
            if not paths[0].exists():
                print(f"  SKIP {tissue}: {paths[0].name} not found.")
            else:
                tissue_bedgraphs[tissue] = paths[0]
                print(f"  OK   {tissue}: {paths[0].name}")
        else:
            merged = BGDIR / f"{tissue}_merged_RC.bedgraph"
            missing = [p for p in paths if not p.exists()]
            if missing:
                print(f"  SKIP {tissue}: missing {[p.name for p in missing]}")
            else:
                print(f"  {tissue}:", end=" ")
                if merge_bedgraphs(paths, merged):
                    tissue_bedgraphs[tissue] = merged

    available = sorted(tissue_bedgraphs)
    n_tissues = len(available)
    print(f"\n  {n_tissues} tissues available: {available}")

    if n_tissues == 0:
        print("ERROR: no tissue bedgraphs found.")
        windows_bed.unlink()
        return

    # ── Step 3: Count signal per window per tissue (parallel) ─────────────────
    print(f"\n[Step 3] Summing bedgraph signal per window "
          f"({n_tissues} tissues, {args.cpus} workers) …")

    raw_cols: dict[str, np.ndarray] = {}

    n_windows = len(windows_filt)

    def _count_tissue(tissue: str) -> tuple[str, np.ndarray]:
        print(f"    {tissue} …")
        return tissue, bedtools_map(windows_bed, tissue_bedgraphs[tissue], n_windows)

    with ThreadPoolExecutor(max_workers=args.cpus) as pool:
        futures = {pool.submit(_count_tissue, t): t for t in available}
        for fut in as_completed(futures):
            tissue, counts = fut.result()
            raw_cols[tissue] = counts
            print(f"    done: {tissue}  "
                  f"(max={counts.max():.4f}, median={np.median(counts):.4f})")

    windows_bed.unlink()

    # Build (n_windows, n_tissues) raw matrix in deterministic column order
    raw_matrix = np.column_stack([raw_cols[t] for t in available])  # (W, T)

    # ── Step 4: QQ normalisation ───────────────────────────────────────────────
    print("\n[Step 4] QQ normalisation …")
    norm_matrix, scale_factors = qq_normalize(raw_matrix, floor=args.floor)

    print("  Scale factors:")
    for tissue, sf in zip(available, scale_factors):
        print(f"    {tissue:25s}  {sf:.6f}")

    # Save the normalised count matrix for inspection
    norm_df = pd.DataFrame(
        norm_matrix,
        columns=available,
    )
    norm_df.insert(0, "type", windows_filt["type"].values)
    norm_df.insert(0, "end",  windows_filt["end"].values)
    norm_df.insert(0, "start",windows_filt["start"].values)
    norm_df.insert(0, "chr",  windows_filt["chr"].values)
    norm_tsv = OUT_DIR / "normalised_counts.tsv.gz"
    norm_df.to_csv(norm_tsv, sep="\t", index=False, compression="gzip")
    print(f"  Normalised counts → {norm_tsv}")

    # ── Step 5: Marker selection ──────────────────────────────────────────────
    top_n = args.top_n
    print(f"\n[Step 5] Marker selection "
          f"(focal > {args.high}, all others < {args.low}"
          + (f", top-{top_n} by SNR" if top_n > 0 else "") + ") …")

    tissue_idx  = {t: i for i, t in enumerate(available)}
    summary_rows = []
    snr_curves   = {}   # tissue → sorted SNR array of all passing markers (for plot)

    for tissue in available:
        j      = tissue_idx[tissue]
        focal  = norm_matrix[:, j]
        others = np.delete(norm_matrix, j, axis=1)   # (W, T-1)
        other_max = others.max(axis=1)

        is_marker = (focal > args.high) & (other_max < args.low)
        pass_idx  = np.where(is_marker)[0]

        # SNR = focal / max(others) for every passing window
        snr = focal[pass_idx] / (other_max[pass_idx] + 1e-6)
        snr_curves[tissue] = np.sort(snr)[::-1]        # descending, for plot

        # Apply top-N trim if requested
        if top_n > 0 and len(pass_idx) > top_n:
            keep = np.argsort(snr)[::-1][:top_n]
            pass_idx = pass_idx[keep]
            snr      = snr[keep]

        marker_df = windows_filt.iloc[pass_idx].copy()
        marker_df["focal_score"] = focal[pass_idx]
        marker_df["snr"]         = snr

        n_total = len(marker_df)
        n_tss   = (marker_df["type"] == "TSS").sum()
        n_extra = (marker_df["type"] == "EXTRA_TSS").sum()
        n_enh   = (marker_df["type"] == "Enhancer").sum()

        bed_path = OUT_DIR / f"{tissue}_markers.bed.gz"
        with gzip.open(bed_path, "wt") as fh:
            if n_total > 0:
                out = marker_df[["chr", "start", "end", "type"]].copy()
                out["score"] = marker_df["focal_score"].round(4)
                out["snr"]   = marker_df["snr"].round(4)
                fh.write(out.to_csv(sep="\t", header=False, index=False))

        n_pass_before_trim = len(snr_curves[tissue])
        summary_rows.append({
            "tissue":         tissue,
            "n_pass_thresh":  n_pass_before_trim,
            "n_markers":      n_total,
            "n_TSS":          n_tss,
            "n_EXTRA_TSS":    n_extra,
            "n_Enhancer":     n_enh,
            "snr_median":     float(np.median(snr)) if n_total else 0.0,
            "bed":            str(bed_path),
        })

        trim_note = (f"  (trimmed from {n_pass_before_trim:,})"
                     if top_n > 0 and n_pass_before_trim > top_n else "")
        print(f"  {tissue:25s}  {n_total:6,} markers{trim_note}  "
              f"(TSS {n_tss:,}  EXTRA_TSS {n_extra:,}  Enhancer {n_enh:,})")

    # ── Step 5b: SNR rank-curve diagnostics plot ──────────────────────────────
    plot_path = OUT_DIR / "snr_rank_curves.png"
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.get_cmap("tab20", len(snr_curves))
    for idx, (tissue, curve) in enumerate(sorted(snr_curves.items())):
        ranks = np.arange(1, len(curve) + 1)
        ax.plot(ranks, curve, lw=1.2, label=tissue, color=cmap(idx))
    if top_n > 0:
        ax.axvline(top_n, color="black", ls="--", lw=1.2, label=f"top-N = {top_n}")
    ax.set_xlabel("Rank (by SNR, descending)")
    ax.set_xlim(xmin=0, xmax=2000)
    ax.set_ylabel("SNR  [focal / max(others)]")
    ax.set_title("SNR rank curves per tissue — all markers passing hard thresholds")
    ax.set_yscale("log")
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  SNR rank curves → {plot_path}")

    # ── Step 6: Summary ───────────────────────────────────────────────────────
    summary_df = pd.DataFrame(summary_rows,
                              columns=["tissue", "n_pass_thresh", "n_markers",
                                       "n_TSS", "n_EXTRA_TSS", "n_Enhancer",
                                       "snr_median", "bed"])
    summary_tsv = OUT_DIR / "marker_summary.tsv"
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    print(f"\n[Summary] → {summary_tsv}")
    print(summary_df.to_string(index=False))
    print(f"\nDone. Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
