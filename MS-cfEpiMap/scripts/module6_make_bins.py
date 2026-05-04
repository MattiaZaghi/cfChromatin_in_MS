#!/usr/bin/env python3
"""
module6_make_bins.py — Generate genome-tiling 5 kb bins.

PURPOSE
Tile the entire hg19 genome into fixed-width bins to complement the
reference-guided RRE analysis with an unbiased, annotation-independent view.
A differential bin can flag regulatory activity differences that are not
captured by any published H3K27ac reference peak set (e.g. patient-specific
neo-enhancers, novel regulatory elements in cfDNA).

After tiling:
  1. Remove bins overlapping DAC exclusion regions (mappability artefacts).
  2. Optionally restrict to bins that overlap the RRE universe — this reduces
     the multiple-testing burden (~10× fewer features) while retaining all
     biologically annotated regions. Recommended when sample sizes are small.

The output BED is used in Module 6 (bin_counts.py) as the counting interval
file, and later by DESeq2 (Module 7c) for bin-level differential analysis.

Usage (called by Snakemake module6_bins.smk):
    python scripts/module6_make_bins.py \\
        --chrom_sizes reference/genome/hg19.chrom.sizes \\
        --dac reference/genome/hg19_dac_exclusion.bed \\
        --rre reference/rre/ms_rre_universe.bed \\
        --bin_size 5000 \\
        --restrict_rre true \\
        --output results/counts/bin_universe.bed
"""

import argparse
import subprocess
import tempfile
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom_sizes",  required=True,
                    help="hg19.chrom.sizes file (two columns: chr  length)")
    ap.add_argument("--dac",          required=True,
                    help="ENCODE DAC exclusion BED — removes mappability artefacts")
    ap.add_argument("--rre",          required=True,
                    help="RRE universe BED — used for optional overlap restriction")
    ap.add_argument("--bin_size",     type=int, default=5000,
                    help="Bin width in bp (default 5000)")
    ap.add_argument("--restrict_rre", default="true",
                    help="'true' to keep only bins that overlap the RRE universe")
    ap.add_argument("--output",       required=True,
                    help="Path for the output bin BED file")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # Temporary files for intermediate steps, cleaned up in the finally block
    with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as tmp_bins, \
         tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as tmp_nodac:
        bins_path  = tmp_bins.name
        nodac_path = tmp_nodac.name

    try:
        # ── Step 1: Tile the genome into fixed-width windows ─────────────────
        # bedtools makewindows divides every chromosome into non-overlapping
        # bins of exactly bin_size bp (the last bin on each chromosome will
        # be shorter). Output: 3-column BED (chr, start, end).
        cmd = f"bedtools makewindows -g {args.chrom_sizes} -w {args.bin_size}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"makewindows error: {result.stderr}", file=sys.stderr)
            sys.exit(1)
        with open(bins_path, "w") as f:
            f.write(result.stdout)

        # ── Step 2: Remove DAC exclusion regions ─────────────────────────────
        # bedtools subtract -A removes any bin that has even 1 bp overlap with
        # a DAC region. DAC regions are ENCODE's "Difficult-to-map Accessible
        # Chromatin" list: repetitive, low-complexity, or blacklisted loci where
        # short-read alignment is unreliable, producing spurious read pileups.
        if os.path.exists(args.dac):
            cmd = f"bedtools subtract -a {bins_path} -b {args.dac} -A"
        else:
            print("[Module 6] WARNING: DAC file not found, skipping DAC removal.",
                  file=sys.stderr)
            cmd = f"cat {bins_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        with open(nodac_path, "w") as f:
            f.write(result.stdout)

        # ── Step 3: Optionally restrict to RRE-overlapping bins ──────────────
        # bedtools intersect -u keeps only bins that have at least 1 bp overlap
        # with the RRE universe. This is a trade-off:
        #   Pro: fewer tests → higher statistical power per test, faster runtime.
        #   Con: misses regulatory elements absent from the reference panel.
        # Recommended for small cohorts (< 20 samples per group).
        if str(args.restrict_rre).lower() in ("true", "1", "yes"):
            cmd = f"bedtools intersect -a {nodac_path} -b {args.rre} -u"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            final = result.stdout
        else:
            # No restriction: use all DAC-filtered bins (genome-wide unbiased scan)
            with open(nodac_path) as f:
                final = f.read()

        # ── Step 4: Write the final bin BED ──────────────────────────────────
        with open(args.output, "w") as f:
            f.write(final)

        n = final.count("\n")
        print(f"[Module 6] {n} bins written to {args.output}")

    finally:
        # Clean up temp files regardless of success or failure
        for p in (bins_path, nodac_path):
            if os.path.exists(p):
                os.remove(p)


if __name__ == "__main__":
    main()
