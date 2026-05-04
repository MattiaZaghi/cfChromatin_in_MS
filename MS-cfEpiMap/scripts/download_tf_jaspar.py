#!/usr/bin/env python3
"""
Download JASPAR 2022 CORE PWMs and scan for TF binding sites.
No FIMO / MEME suite required — uses a numpy PWM scanner.

Requirements (pip or conda):
  requests  numpy  (bedtools must be in PATH for --region-bed mode)

Usage (from MS-cfEpiMap/)
--------------------------
  # Scan only RRE universe (~5-15 min total, recommended):
  python scripts/download_tf_jaspar.py \
      --region-bed reference/rre/ms_rre_universe.bed

  # Scan full hg19 (slow, ~30 min–1 h total):
  python scripts/download_tf_jaspar.py
"""

import argparse
import io
import os
import re
import subprocess
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import requests

# ── TFs to process ────────────────────────────────────────────────────────────
TF_NAMES = [
    "PAX5",  "IRF4",  "EBF1",          # B cell
    "TBX21", "RORC",  "FOXP3",         # T cell / MS inflammation
    "OLIG2", "MYRF",                   # CNS / oligodendrocyte
    "SPI1",  "IRF8",                   # myeloid / innate
    "STAT1", "IRF1",  "RELA",  "YY1",  # broad reference
]

# ── Paths ─────────────────────────────────────────────────────────────────────
GENOME_FA  = "/home/mattia/Genomes/hg19/fa/hg19.fa"
OUTPUT_DIR = Path("reference/tf_sites")
PWM_DIR    = Path("reference/tf_pwms")

# JASPAR 2022 CORE vertebrates – single ZIP (~5 MB)
JASPAR_ZIP_URL = (
    "https://jaspar.elixir.no/download/data/2022/CORE/"
    "JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.zip"
)

# PWM score threshold: fraction of maximum possible score
# 0.80 = fairly strict; lower it to 0.75 if too few sites
SCORE_FRACTION = 0.80
MIN_SITES      = 100


# ── JASPAR parsing ────────────────────────────────────────────────────────────

def download_all_pfms() -> dict[str, list[dict]]:
    """Download JASPAR 2022 ZIP and parse every vertebrate PFM.

    Returns {TF_NAME_UPPER: [pfm_dict, ...]} where pfm_dict has keys A C G T.
    """
    print("Downloading JASPAR 2022 CORE vertebrates PFMs ...")
    resp = requests.get(JASPAR_ZIP_URL, timeout=120)
    resp.raise_for_status()
    print(f"  Downloaded {len(resp.content) / 1e6:.1f} MB")

    matrices: dict[str, list[dict]] = {}
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        for fname in zf.namelist():
            if not fname.endswith(".jaspar"):
                continue
            tf_name, pfm = _parse_jaspar_file(zf.read(fname).decode())
            if tf_name:
                matrices.setdefault(tf_name.upper(), []).append(pfm)

    print(f"  Parsed {sum(len(v) for v in matrices.values())} matrices "
          f"for {len(matrices)} TFs")
    return matrices


def _parse_jaspar_file(content: str):
    """Parse a single .jaspar PFM block. Returns (tf_name, pfm_dict)."""
    lines = content.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        return None, None

    # Header: ">MA0014.3 PAX5"
    parts = lines[0][1:].strip().split()
    tf_name = parts[1] if len(parts) > 1 else parts[0]

    pfm: dict[str, np.ndarray] = {}
    for line in lines[1:]:
        m = re.match(r"([ACGT])\s*\[?\s*([\d\s.]+?)\s*\]?\s*$", line.strip())
        if m:
            pfm[m.group(1)] = np.array(m.group(2).split(), dtype=float)

    return (tf_name, pfm) if len(pfm) == 4 else (None, None)


def pfm_to_pwm(pfm: dict) -> np.ndarray:
    """Convert PFM (counts) to log2 PWM (log-odds vs uniform background).

    Returns ndarray of shape (4, motif_length), rows ordered A C G T.
    """
    mat = np.array([pfm[n] for n in "ACGT"], dtype=float)  # (4, L)
    mat += 0.1          # Laplace pseudocount
    mat /= mat.sum(axis=0, keepdims=True)   # PPM
    mat = np.log2(mat / 0.25)              # log-odds vs uniform background
    return mat


def best_pwm(pfm_list: list[dict]) -> np.ndarray:
    """Pick the matrix with the highest information content."""
    pwms = [pfm_to_pwm(p) for p in pfm_list]
    ic   = [pwm.max(axis=0).sum() for pwm in pwms]
    return pwms[int(np.argmax(ic))]


# ── Sequence scanner ──────────────────────────────────────────────────────────

_NUC_MAP = np.full(256, -1, dtype=np.int8)
for _i, _c in enumerate(b"ACGTacgt"):
    _NUC_MAP[_c] = _i % 4


def _scan_single(seq: str, pwm: np.ndarray, threshold: float) -> np.ndarray:
    """Return 0-based start positions where PWM score >= threshold."""
    L = pwm.shape[1]
    n = len(seq)
    if n < L:
        return np.array([], dtype=np.int64)

    seq_idx = _NUC_MAP[np.frombuffer(seq.encode(), dtype=np.uint8)]
    scores  = np.zeros(n - L + 1, dtype=np.float32)

    for pos in range(L):
        sl = seq_idx[pos : n - L + pos + 1]
        ok = sl >= 0
        if ok.any():
            scores[ok] += pwm[sl[ok], pos]

    return np.where(scores >= threshold)[0].astype(np.int64)


def scan_fasta(fasta_path: str, pwm: np.ndarray, threshold: float,
               used_region_bed: bool) -> list[tuple[str, int, int]]:
    """Scan all sequences in a FASTA file.

    If used_region_bed=True, sequence headers are 'chr:start-end' (bedtools
    getfasta format) and positions are converted to genome coordinates.
    """
    hits: list[tuple[str, int, int]] = []
    L = pwm.shape[1]

    chrom_id  = None
    seq_parts: list[str] = []

    def flush():
        if chrom_id is None:
            return
        seq = "".join(seq_parts).upper()
        positions = _scan_single(seq, pwm, threshold)

        if used_region_bed and ":" in chrom_id:
            chrom, coords = chrom_id.rsplit(":", 1)
            offset = int(coords.split("-")[0])
        else:
            chrom  = chrom_id
            offset = 0

        for p in positions:
            hits.append((chrom, offset + int(p), offset + int(p) + L))

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                flush()
                chrom_id  = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        flush()

    return hits


# ── Helpers ───────────────────────────────────────────────────────────────────

def extract_fasta(region_bed: str, genome_fa: str, tmp_path: str) -> None:
    subprocess.run(
        ["bedtools", "getfasta", "-fi", genome_fa, "-bed", region_bed,
         "-fo", tmp_path],
        check=True,
    )


def write_bed(hits: list[tuple[str, int, int]], out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for chrom, start, end in hits:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    return len(hits)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate TF binding site BEDs from JASPAR 2022 (no FIMO needed)."
    )
    parser.add_argument(
        "--region-bed", default=None, metavar="BED",
        help="Scan only these regions (recommended: reference/rre/ms_rre_universe.bed). "
             "Much faster than full genome.",
    )
    parser.add_argument(
        "--genome",   default=GENOME_FA,
        help=f"hg19 FASTA  [default: {GENOME_FA}]",
    )
    parser.add_argument(
        "--tfs", nargs="+", default=TF_NAMES,
        help="TF names to process (default: all priority TFs)",
    )
    parser.add_argument(
        "--score-fraction", type=float, default=SCORE_FRACTION,
        help=f"Min PWM score as fraction of max possible score  [default: {SCORE_FRACTION}]",
    )
    parser.add_argument(
        "--output-dir", default=str(OUTPUT_DIR),
        help=f"Where to write TF BED files  [default: {OUTPUT_DIR}]",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    PWM_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Download all PFMs ──────────────────────────────────────────────────
    all_pfms = download_all_pfms()

    # ── 2. Prepare target FASTA ───────────────────────────────────────────────
    tmp_fa   = None
    if args.region_bed:
        tmp_fa = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        print(f"\nExtracting sequences for {args.region_bed} ...")
        extract_fasta(args.region_bed, args.genome, tmp_fa)
        print(f"  Saved to {tmp_fa}")
        target_fa = tmp_fa
    else:
        target_fa = args.genome

    # ── 3. Scan each TF ───────────────────────────────────────────────────────
    print(f"\nScanning {len(args.tfs)} TFs "
          f"(score >= {args.score_fraction:.0%} of max) ...\n")
    skipped = []

    for tf in args.tfs:
        tf_key = tf.upper()
        if tf_key not in all_pfms:
            print(f"[{tf}]  SKIP — not found in JASPAR 2022 CORE vertebrates")
            skipped.append(tf)
            continue

        pfm_list = all_pfms[tf_key]
        pwm      = best_pwm(pfm_list)
        max_score = pwm.max(axis=0).sum()
        threshold = args.score_fraction * max_score

        # Save PWM as text for reference
        np.savetxt(PWM_DIR / f"{tf}.pwm.txt", pwm,
                   header=f"{tf}  threshold={threshold:.3f}", fmt="%.4f")

        print(f"[{tf}]  motif_len={pwm.shape[1]}  "
              f"max_score={max_score:.2f}  threshold={threshold:.2f}")

        hits   = scan_fasta(target_fa, pwm, threshold, args.region_bed is not None)
        n_hits = write_bed(hits, out_dir / f"{tf}.bed")

        if n_hits < MIN_SITES:
            print(f"  WARNING: {n_hits} sites — below MIN_SITES={MIN_SITES}. "
                  f"Try --score-fraction 0.75")
        else:
            print(f"  {n_hits:,} sites  ->  {out_dir}/{tf}.bed")

    # ── 4. Cleanup ────────────────────────────────────────────────────────────
    if tmp_fa:
        os.remove(tmp_fa)

    print(f"\nDone. BED files in {out_dir}/")
    if skipped:
        print(f"Skipped (not in JASPAR 2022 CORE vertebrates): {skipped}")


if __name__ == "__main__":
    main()
