#!/usr/bin/env python3
"""
Build H3K27ac on-target and off-target BED files from Roadmap Epigenomics
15-state chromHMM annotations.

On-target  : 200 bp windows annotated as H3K27ac-active (1_TssA, 2_TssAFlnk,
             6_EnhG, 7_Enh) in >= 95% of samples.
             After merging: keep regions > 1 kb, remove DAC/blacklist, chrX/Y.

Off-target : 200 bp windows NEVER annotated as active in any sample.
             Remove windows within 10 kb of any on-target site.
             After merging: keep regions > 1 kb, remove DAC/blacklist, chrX/Y.

Usage:
    python build_chromhmm_reference.py \
        --chromhmm_dir  /date/gcb/gcb_MZ/Chrom_HMM_hg19 \
        --chrom_sizes   /home/mattia/Genomes/hg19/hg19.chrom.sizes \
        --blacklist     /home/mattia/Genomes/hg19/hg19-blacklist.v2.bed \
        --out_dir       ref_files \
        --threads       8
"""

import argparse
import glob
import os
import subprocess
import sys
import tempfile

ACTIVE_STATES       = {"1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh"}
MIN_REGION_BP       = 1000
OFF_TARGET_FLANK_KB = 10000
AUTOSOME_CHROMS     = {f"chr{i}" for i in range(1, 23)}


# ── helpers ───────────────────────────────────────────────────────────────────

def run(cmd, desc=""):
    print(f"  [RUN] {desc or cmd[:80]}", flush=True)
    result = subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
    if result.stderr:
        msg = result.stderr.decode().strip()
        if msg:
            print(f"  [STDERR] {msg[:300]}", flush=True)


def count_lines(path):
    r = subprocess.run(f"wc -l < {path}", shell=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return int(r.stdout.decode().strip() or 0)


def sum_bp(path):
    r = subprocess.run(f"awk '{{s+=$3-$2}}END{{print s+0}}' {path}", shell=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return int(r.stdout.decode().strip() or 0)


def count_chromhmm_files(chromhmm_dir):
    files = sorted(glob.glob(
        os.path.join(chromhmm_dir, "E*_15_coreMarks_mnemonics.bed")))
    if not files:
        sys.exit(f"ERROR: no chromHMM BED files found in {chromhmm_dir}")
    return files


# ── Step 1 ────────────────────────────────────────────────────────────────────

def build_active_sorted_bed(files, tmpdir, threads):
    """
    Read all chromHMM files in Python (no shell quoting issues),
    filter for active states on autosomes, write a sorted BED.
    Returns path to sorted BED and sample count.
    """
    n = len(files)
    print(f"  Reading {n} chromHMM files (Python, no shell quoting issues)...",
          flush=True)

    raw_bed = os.path.join(tmpdir, "all_active_raw.bed")
    written = 0
    with open(raw_bed, "w") as out:
        for f in files:
            with open(f) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    chrom, start, end, state = parts[0], parts[1], parts[2], parts[3]
                    if chrom in AUTOSOME_CHROMS and state in ACTIVE_STATES:
                        out.write(f"{chrom}\t{start}\t{end}\n")
                        written += 1

    print(f"  Active rows written: {written:,}", flush=True)

    sorted_bed = os.path.join(tmpdir, "all_active_sorted.bed")
    run(f"sort --parallel={threads} -k1,1 -k2,2n {raw_bed} > {sorted_bed}",
        "Sorting active regions")
    print(f"  sorted.bed lines: {count_lines(sorted_bed):,}", flush=True)
    return sorted_bed, n


# ── Step 2 ────────────────────────────────────────────────────────────────────

def build_coverage(sorted_bed, chrom_sizes, tmpdir):
    coverage_bg = os.path.join(tmpdir, "active_coverage.bg")
    run(f"bedtools genomecov -i {sorted_bed} -g {chrom_sizes} -bg > {coverage_bg}",
        "Computing per-position activity depth")
    print(f"  coverage.bg lines: {count_lines(coverage_bg):,}", flush=True)
    return coverage_bg


# ── Step 3 ────────────────────────────────────────────────────────────────────

def build_ontarget(coverage_bg, n, tmpdir, blacklist, chrom_sizes, out_dir,
                   fraction_threshold=0.0):
    threshold = int(n * fraction_threshold) + 1
    print(f"  On-target threshold: >= {threshold}/{n} samples "
          f"({fraction_threshold*100:.0f}%)", flush=True)

    high_conf = os.path.join(tmpdir, "high_conf_active.bed")
    run(f"awk '$4 >= {threshold}' {coverage_bg} | cut -f1-3 > {high_conf}",
        f"Filtering depth >= {threshold}")
    print(f"  high_conf.bed: {count_lines(high_conf):,} lines", flush=True)

    merged = os.path.join(tmpdir, "ontarget_merged.bed")
    run(f"bedtools merge -i {high_conf} > {merged}", "Merging on-target")
    print(f"  merged.bed: {count_lines(merged):,} lines", flush=True)

    sized = os.path.join(tmpdir, "ontarget_sized.bed")
    run(f"awk '$3-$2 >= {MIN_REGION_BP}' {merged} > {sized}", "Size filter on-target")
    print(f"  sized.bed: {count_lines(sized):,} lines", flush=True)

    out = os.path.join(out_dir, "chromhmm_ontarget_H3K27ac.bed")
    if blacklist and os.path.getsize(blacklist) > 0:
        run(f"bedtools intersect -a {sized} -b {blacklist} -v | "
            f"bedtools sort -g {chrom_sizes} -i - > {out}",
            "Removing blacklist, sorting on-target")
    else:
        run(f"bedtools sort -g {chrom_sizes} -i {sized} > {out}",
            "Sorting on-target")

    n_regions = count_lines(out)
    bp_total  = sum_bp(out)
    print(f"  On-target: {n_regions:,} regions, {bp_total:,} bp", flush=True)
    return out


# ── Step 4 ────────────────────────────────────────────────────────────────────

def build_offtarget(sorted_active_bed, ontarget_bed, chrom_sizes,
                    blacklist, tmpdir, out_dir):
    # Autosome-only chrom sizes (written first, used for sorting)
    auto_sizes = os.path.join(tmpdir, "autosomes.sizes")
    with open(chrom_sizes) as fh, open(auto_sizes, "w") as out:
        for line in fh:
            chrom = line.split("\t")[0]
            if chrom in AUTOSOME_CHROMS:
                out.write(line)

    # Ever-active union — sort against genome file so complement works
    ever_active = os.path.join(tmpdir, "ever_active_merged.bed")
    run(f"bedtools merge -i {sorted_active_bed} | "
        f"bedtools sort -g {auto_sizes} -i - > {ever_active}",
        "Building ever-active union (genome-sorted)")

    # Complement = never-active
    never_active = os.path.join(tmpdir, "never_active.bed")
    run(f"bedtools complement -i {ever_active} -g {auto_sizes} > {never_active}",
        "Building never-active complement")
    print(f"  never_active.bed: {count_lines(never_active):,} lines", flush=True)

    # Expand on-target by OFF_TARGET_FLANK_KB
    ontarget_flanked = os.path.join(tmpdir, "ontarget_flanked.bed")
    run(f"bedtools slop -i {ontarget_bed} -g {chrom_sizes} "
        f"-b {OFF_TARGET_FLANK_KB} | bedtools merge -i - > {ontarget_flanked}",
        f"Expanding on-target by {OFF_TARGET_FLANK_KB//1000} kb")

    no_flank = os.path.join(tmpdir, "offtarget_no_flank.bed")
    run(f"bedtools intersect -a {never_active} -b {ontarget_flanked} -v "
        f"> {no_flank}",
        "Removing on-target flanks")

    merged = os.path.join(tmpdir, "offtarget_merged.bed")
    run(f"bedtools merge -i {no_flank} > {merged}", "Merging off-target")

    sized = os.path.join(tmpdir, "offtarget_sized.bed")
    run(f"awk '$3-$2 >= {MIN_REGION_BP}' {merged} > {sized}",
        "Size filter off-target")
    print(f"  sized.bed: {count_lines(sized):,} lines", flush=True)

    out = os.path.join(out_dir, "chromhmm_offtarget_H3K27ac.bed")
    if blacklist and os.path.getsize(blacklist) > 0:
        run(f"bedtools intersect -a {sized} -b {blacklist} -v | "
            f"bedtools sort -g {chrom_sizes} -i - > {out}",
            "Removing blacklist, sorting off-target")
    else:
        run(f"bedtools sort -g {chrom_sizes} -i {sized} > {out}",
            "Sorting off-target")

    n_regions = count_lines(out)
    bp_total  = sum_bp(out)
    print(f"  Off-target: {n_regions:,} regions, {bp_total:,} bp", flush=True)
    return out


# ── metadata ──────────────────────────────────────────────────────────────────

def write_metadata(out_dir, n_files, ontarget, offtarget, chrom_sizes, blacklist,
                   fraction_threshold=0.0):
    meta = os.path.join(out_dir, "chromhmm_reference_metadata.txt")
    with open(meta, "w") as fh:
        fh.write(f"mark=H3K27ac\n")
        fh.write(f"chromhmm_model=15-state Roadmap Epigenomics core marks\n")
        fh.write(f"n_samples={n_files}\n")
        fh.write(f"active_states={','.join(sorted(ACTIVE_STATES))}\n")
        fh.write(f"ontarget_fraction_threshold={fraction_threshold}\n")
        fh.write(f"offtarget_flank_exclusion_bp={OFF_TARGET_FLANK_KB}\n")
        fh.write(f"min_region_bp={MIN_REGION_BP}\n")
        fh.write(f"blacklist={blacklist}\n")
        fh.write(f"chrom_sizes={chrom_sizes}\n")
        fh.write(f"ontarget_n_regions={count_lines(ontarget)}\n")
        fh.write(f"ontarget_total_bp={sum_bp(ontarget)}\n")
        fh.write(f"offtarget_n_regions={count_lines(offtarget)}\n")
        fh.write(f"offtarget_total_bp={sum_bp(offtarget)}\n")
    print(f"  Metadata written to {meta}", flush=True)


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--chromhmm_dir", required=True)
    parser.add_argument("--chrom_sizes",  required=True)
    parser.add_argument("--blacklist",    default="")
    parser.add_argument("--out_dir",      default="ref_files")
    parser.add_argument("--threads",      type=int, default=8)
    parser.add_argument("--tmpdir",       default=None)
    parser.add_argument(
        "--ontarget_threshold", type=float, default=0.0,
        help=(
            "Fraction of samples that must have a region annotated as active "
            "for it to be included in the on-target set.  "
            "0.0 (default) = union: active in >=1 sample (~300-400k regions). "
            "0.95 = active in >=95%% of samples (~4k constitutively-active loci)."
        )
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    files = count_chromhmm_files(args.chromhmm_dir)
    print(f"Found {len(files)} chromHMM files", flush=True)

    with tempfile.TemporaryDirectory(
            dir=args.tmpdir or "/tmp", prefix="chromhmm_ref_") as tmpdir:

        print("\n=== Step 1: Extracting active regions ===", flush=True)
        sorted_active_bed, n = build_active_sorted_bed(files, tmpdir, args.threads)

        print("\n=== Step 2: Computing coverage depth ===", flush=True)
        coverage_bg = build_coverage(sorted_active_bed, args.chrom_sizes, tmpdir)

        print("\n=== Step 3: Building on-target BED ===", flush=True)
        ontarget = build_ontarget(
            coverage_bg, n, tmpdir, args.blacklist, args.chrom_sizes, args.out_dir,
            fraction_threshold=args.ontarget_threshold)

        print("\n=== Step 4: Building off-target BED ===", flush=True)
        offtarget = build_offtarget(
            sorted_active_bed, ontarget, args.chrom_sizes,
            args.blacklist, tmpdir, args.out_dir)

        write_metadata(args.out_dir, n, ontarget, offtarget,
                       args.chrom_sizes, args.blacklist,
                       fraction_threshold=args.ontarget_threshold)

    print("\nDone.", flush=True)
