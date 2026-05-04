#!/usr/bin/env python3
"""
count_midpoints.py — Count fragment midpoints overlapping genomic regions.

Usage:
    python scripts/count_midpoints.py \
        --frags   sample.bed \
        --regions regions.bed \
        --output  counts.bed

Output: same columns as regions.bed plus a final integer count column.
Called by multiple modules as a core counting utility.

WHY MIDPOINTS?
Each row in the fragment BED represents one complete DNA fragment (chr, start, end).
The fragment midpoint (start + end) // 2 is used as the single representative
position because:
  1. It best approximates the nucleosomal dyad under histone-bound cfDNA,
     where H3K27ac enrichment peaks at the centre of the marked element.
  2. It avoids double-counting: a long fragment bridging a narrow region boundary
     would be counted at most once, from wherever its midpoint falls.
  3. It is the approach used in the Sadeh et al. cfChIP-seq framework on which
     this pipeline is based.
"""

import argparse
import subprocess
import tempfile
import os
import sys


def count_midpoints(frags_bed: str, regions_bed: str, output_file: str) -> None:
    """
    Count fragment midpoints overlapping each region in regions_bed.

    Steps:
      1. Convert every fragment (chr, start, end) → 1-bp midpoint interval
         (chr, mid, mid+1) using awk.
      2. Sort the midpoints BED so bedtools intersect can use its sweep-line
         algorithm (required for correct counting).
      3. Run bedtools intersect -c to count how many midpoints fall inside
         each region; the count is appended as the last column.
      4. Write the result to output_file.

    A temporary file is used for the sorted midpoints so we do not need to
    hold everything in memory and the file is cleaned up on exit.
    """

    # Create a temp file in the same directory as the output to avoid
    # cross-device rename issues on network filesystems.
    with tempfile.NamedTemporaryFile(
        suffix=".bed", delete=False, mode="w", dir=os.path.dirname(output_file) or "."
    ) as tmp:
        tmp_path = tmp.name

    try:
        # ── Step 1 + 2: compute midpoints and sort in a single pipeline ──────
        # awk computes mid = int((start+end)/2) and emits a 1-bp interval.
        # bedtools sort -i - reads from stdin and sorts by chromosome then position.
        # The sorted output is written to tmp_path for use in the intersect step.
        awk = (
            "awk 'BEGIN{OFS=\"\\t\"} "
            "{mid=int(($2+$3)/2); print $1, mid, mid+1}' "
            f"{frags_bed}"
        )
        sort = "bedtools sort -i -"
        with open(tmp_path, "w") as out:
            p1 = subprocess.Popen(awk, shell=True, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                sort, shell=True, stdin=p1.stdout, stdout=out, stderr=subprocess.PIPE
            )
            p1.stdout.close()   # allow p1 to receive SIGPIPE if p2 exits early
            _, err = p2.communicate()
            if p2.returncode != 0:
                print(f"bedtools sort error: {err.decode()}", file=sys.stderr)
                sys.exit(1)

        # ── Step 3: count midpoints per region with bedtools intersect -c ────
        # -a: the query regions (each line gets a count appended)
        # -b: the midpoints we just computed
        # -c: count the number of -b entries that overlap each -a entry
        cmd = f"bedtools intersect -a {regions_bed} -b {tmp_path} -c"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"bedtools intersect error: {result.stderr}", file=sys.stderr)
            sys.exit(1)

        # ── Step 4: write counts to output file ───────────────────────────────
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        with open(output_file, "w") as fh:
            fh.write(result.stdout)

    finally:
        # Always clean up the temporary sorted-midpoints file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count fragment midpoints in regions.")
    parser.add_argument("--frags",   required=True, help="Fragment BED (3-column)")
    parser.add_argument("--regions", required=True, help="Regions BED file")
    parser.add_argument("--output",  required=True, help="Output file path")
    args = parser.parse_args()
    count_midpoints(args.frags, args.regions, args.output)
