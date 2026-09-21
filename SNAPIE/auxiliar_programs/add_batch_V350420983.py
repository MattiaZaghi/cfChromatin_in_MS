#!/usr/bin/env python3
"""
Register the V350420983 fastq batch: create symlinks in the shared fastq
directory and append rows to the SNAPIE samplesheet.

Group assignment for this batch (stage/treatment largely unknown):
  H*      -> Ctrl                (healthy donors)
  Nat-*   -> MS-Natalizumab      (natalizumab-treated MS)
  <digits>-> MS-Unknown          (coded samples, assumed MS)

Run with --dry-run to print the mapping without writing anything.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path

SRC_DIR = Path("/date/gcb/gcb_MZ/XG_cfChromatin_Frozen/V350420983")
LINK_DIR = Path("/date/gcb/gcb_MZ/cfChrom_fastqs")
SAMPLESHEET = Path(__file__).resolve().parent.parent / "config" / "samplesheet.csv"

SAMPLE_TYPE = "P"      # plasma
VISIT = "V4"           # protocol/visit tag
INPUT_VOL = "1000"     # plasma input (ul)

SKIP = {"Undetermined"}


def sanitise(name: str) -> str:
    """'#' is hostile to Snakemake wildcards and shells; '-2' style replicate
    suffixes are kept but joined with '_' so they do not look like a new field
    in the '-'-delimited sampleId convention."""
    return name.replace("#", "_")


def parse(stem: str, dual_mark=frozenset()):
    """stem -> (sample_id, mark). Returns None for samples to skip.

    The mark normally lives in the enrichment_mark column, not in the
    sampleId. The exception is a donor assayed with both marks in the same
    batch: sampleId is the only key the pipeline has (trim.smk caches rows in
    a dict keyed on it), so an undisambiguated ID would silently drop one of
    the two fastq pairs. Those donors keep a -H3K4me3 suffix.
    """
    if stem in SKIP:
        return None

    mark = "H3K27ac"
    base = stem
    if base.endswith("-H3K4me3"):
        mark = "H3K4me3"
        base = base[: -len("-H3K4me3")]

    if base.startswith("H") and re.fullmatch(r"H\d+", base):
        code, group = base, "Ctrl"
    elif base.lower().startswith("nat-"):
        # Nat-026-ks / Nat-029-KS -> Nat026 (trailing operator initials dropped)
        num = re.match(r"[Nn]at-(\d+)", base).group(1)
        code, group = f"Nat{num}", "MS-Natalizumab"
    else:
        # 03-114#3 -> 03114_3 ; 09-037-2 -> 09037_2 ; 19-104 -> 19104
        m = re.fullmatch(r"(\d+)-(\d+)(?:[#-](\d+))?", sanitise(base).replace("_", "#"))
        if not m:
            raise ValueError(f"unrecognised sample name: {stem}")
        code = m.group(1) + m.group(2)
        if m.group(3):
            code += "_" + m.group(3)
        group = "MS-Unknown"

    sample_id = f"{code}-{SAMPLE_TYPE}-{group}-{VISIT}-{INPUT_VOL}"
    if mark == "H3K4me3" and base in dual_mark:
        sample_id += "-H3K4me3"
    return sample_id, mark


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    stems = sorted({p.name[: -len("_R1.fastq.gz")]
                    for p in SRC_DIR.glob("*_R1.fastq.gz")})

    # donors present in this batch under both marks
    dual_mark = {s[: -len("-H3K4me3")] for s in stems
                 if s.endswith("-H3K4me3") and s[: -len("-H3K4me3")] in stems}

    rows, links, problems = [], [], []
    for stem in stems:
        parsed = parse(stem, dual_mark)
        if parsed is None:
            continue
        sample_id, mark = parsed

        link_names = {}
        for read in ("R1", "R2"):
            src = SRC_DIR / f"{stem}_{read}.fastq.gz"
            if not src.exists():
                problems.append(f"missing mate: {src}")
                continue
            dst = LINK_DIR / f"{sanitise(stem)}_{read}.fastq.gz"
            link_names[read] = dst
            links.append((src, dst))
        if len(link_names) != 2:
            continue

        rows.append({
            "sampleId": sample_id,
            "enrichment_mark": mark,
            "read1": str(link_names["R1"]),
            "read2": str(link_names["R2"]),
            "control": "",
        })

    # ---- collision checks -------------------------------------------------
    existing_ids = []
    with open(SAMPLESHEET, newline="") as fh:
        existing_ids = [r["sampleId"] for r in csv.DictReader(fh)]
    new_ids = [r["sampleId"] for r in rows]

    for sid in sorted(set(new_ids) & set(existing_ids)):
        problems.append(f"sampleId already in samplesheet: {sid}")
    for sid in sorted({s for s in new_ids if new_ids.count(s) > 1}):
        problems.append(f"duplicate sampleId within new batch: {sid}")
    for _, dst in links:
        if dst.exists() and not dst.is_symlink():
            problems.append(f"link target exists and is a real file: {dst}")

    for r in rows:
        print(f"{Path(r['read1']).name[:-len('_R1.fastq.gz')]:<20} -> "
              f"{r['sampleId']:<45} {r['enrichment_mark']}")
    print(f"\n{len(rows)} samples "
          f"({sum(r['enrichment_mark'] == 'H3K27ac' for r in rows)} H3K27ac, "
          f"{sum(r['enrichment_mark'] == 'H3K4me3' for r in rows)} H3K4me3), "
          f"{len(links)} symlinks")

    if problems:
        print("\nPROBLEMS:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print("\n[dry-run] nothing written")
        return

    LINK_DIR.mkdir(parents=True, exist_ok=True)
    for src, dst in links:
        if dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src)
    print(f"\ncreated {len(links)} symlinks in {LINK_DIR}")

    backup = SAMPLESHEET.with_suffix(".csv.bak")
    shutil.copy2(SAMPLESHEET, backup)
    with open(SAMPLESHEET, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sampleId", "enrichment_mark",
                                           "read1", "read2", "control"])
        for r in rows:
            w.writerow(r)
    print(f"appended {len(rows)} rows to {SAMPLESHEET} (backup: {backup.name})")


if __name__ == "__main__":
    main()
