#!/usr/bin/env python3
"""
Dissolve the two duplicate sampleIds in config/samplesheet.csv.

Duplicates are not cosmetic here: trim.smk caches samplesheet rows in a dict
keyed on sampleId, so the second row of a pair silently overwrites the first
and one fastq pair never reaches the pipeline.

1. H17-P-Ctrl-V3-1000 appears twice, once per mark (H3K27ac on P34608_1035,
   H3K4me3 on the hand-delivered H17_K4 run). Same situation as the
   H62/H66/H70/H71/H72 pairs in batch V350420983, so it takes the same
   -H3K4me3 suffix.

2. H6-P-Ctrl-V2-1000 appears twice on two different libraries. The second is
   P34608_1024, which sits between P34608_1023 (H5) and _1025 (H7) -- both of
   which the sheet labels V3. The V2 on that row is a typo; it belongs to the
   V3 series.

NOT touched: rows mapped to P32456_2xxx / P34608_2xxx libraries. Those are
Mattia's own corrections of names mis-submitted to the facility, and the
manifest "User ID" column is the stale side of that comparison. Cross-checking
confirms them -- every 1xxx library they move away from (1031 H14, 1034 H17,
1018 P24136, 1027 18-070) sequenced at 0.00 Mreads.
"""

import csv
import shutil
from pathlib import Path

SHEET = Path(__file__).resolve().parent.parent / "config" / "samplesheet.csv"


def main():
    rows = list(csv.DictReader(open(SHEET)))
    fields = list(rows[0].keys())
    changes = []

    for r in rows:
        sid, mark, r1 = r["sampleId"], r["enrichment_mark"], r["read1"]

        if sid == "H17-P-Ctrl-V3-1000" and mark == "H3K4me3":
            r["sampleId"] = "H17-P-Ctrl-V3-1000-H3K4me3"
        elif sid == "H6-P-Ctrl-V2-1000" and "P34608_1024" in r1:
            r["sampleId"] = "H6-P-Ctrl-V3-1000"
        else:
            continue
        changes.append(f"{sid:<28} -> {r['sampleId']:<30} ({Path(r1).name})")

    for c in changes:
        print("  " + c)

    ids = [r["sampleId"] for r in rows]
    dupes = {s for s in ids if ids.count(s) > 1}
    print(f"\n{len(changes)} rows renamed; {len(rows)} rows, {len(set(ids))} unique sampleIds")
    if dupes:
        raise SystemExit(f"still duplicated: {sorted(dupes)}")

    shutil.copy2(SHEET, SHEET.with_suffix(".csv.predupfix"))
    with open(SHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"written; backup: {SHEET.with_suffix('.csv.predupfix').name}")


if __name__ == "__main__":
    main()
