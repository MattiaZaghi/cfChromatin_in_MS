#!/usr/bin/env python3
"""
Compute per-sample H3K27ac cfChIP-seq enrichment score.

Enrichment score = (on_frags / on_bp) / (off_frags / off_bp)

Fragment midpoints (not read pairs) are intersected against
chromHMM on-target and off-target BED files.
On-target : H3K27ac-active in >= 95% of 131 Roadmap samples
Off-target: never active in any sample, > 10 kb from any on-target site

Usage:
    python compute_enrichment_score.py \
        --frag_bed  results/frags/sample.bed \
        --ontarget  ref_files/chromhmm_ontarget_H3K27ac.bed \
        --offtarget ref_files/chromhmm_offtarget_H3K27ac.bed \
        --sample    sample_id \
        --out       results/enrichment/sample_enrichment.tsv
"""

import argparse
import os
import subprocess
import sys
import tempfile


def run(cmd, check=True):
    result = subprocess.run(cmd, shell=True, check=check,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0 and check:
        sys.exit(f"ERROR running: {cmd}\n{result.stderr.decode()}")
    return result.stdout.decode().strip()


def bed_total_bp(bed_file):
    """Sum of (end - start) across all regions in a BED file."""
    total = run(f"awk '{{s += $3 - $2}} END {{print s+0}}' {bed_file}")
    return int(total)


def count_midpoints_in_bed(frag_bed, target_bed, tmpdir):
    """
    Count fragment midpoints overlapping target_bed.
    Midpoint = int((start + end) / 2)
    """
    midpoints = os.path.join(tmpdir, "midpoints.bed")
    run(f"awk 'OFS=\"\\t\" {{print $1, int(($2+$3)/2), int(($2+$3)/2)+1}}'"
        f" {frag_bed} > {midpoints}")
    count = run(
        f"bedtools intersect -a {midpoints} -b {target_bed} -u | wc -l"
    )
    return int(count)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frag_bed",  required=True,
                        help="3-column fragment BED (chr, start, end)")
    parser.add_argument("--ontarget",  required=True,
                        help="chromHMM on-target BED")
    parser.add_argument("--offtarget", required=True,
                        help="chromHMM off-target BED")
    parser.add_argument("--sample",    required=True, help="Sample ID")
    parser.add_argument("--out",       required=True, help="Output TSV path")
    args = parser.parse_args()

    for f in [args.frag_bed, args.ontarget, args.offtarget]:
        if not os.path.exists(f):
            sys.exit(f"ERROR: file not found: {f}")

    with tempfile.TemporaryDirectory(prefix="enrichment_") as tmpdir:
        print(f"[{args.sample}] Counting midpoints in on-target regions ...", flush=True)
        on_frags = count_midpoints_in_bed(args.frag_bed, args.ontarget, tmpdir)
        on_bp    = bed_total_bp(args.ontarget)

        print(f"[{args.sample}] Counting midpoints in off-target regions ...", flush=True)
        off_frags = count_midpoints_in_bed(args.frag_bed, args.offtarget, tmpdir)
        off_bp    = bed_total_bp(args.offtarget)

    if off_frags == 0 or off_bp == 0:
        print(f"WARNING: [{args.sample}] off-target count is 0; score set to NA", flush=True)
        score = "NA"
    else:
        on_density  = on_frags  / on_bp
        off_density = off_frags / off_bp
        score = round(on_density / off_density, 6)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out, "w") as fh:
        fh.write("sample_id\ton_frags\ton_bp\toff_frags\toff_bp\tenrichment_score\n")
        fh.write(f"{args.sample}\t{on_frags}\t{on_bp}\t{off_frags}\t{off_bp}\t{score}\n")

    print(f"[{args.sample}] Enrichment score = {score}  "
          f"(on: {on_frags:,} frags / {on_bp:,} bp; "
          f"off: {off_frags:,} frags / {off_bp:,} bp)",
          flush=True)


if __name__ == "__main__":
    main()
