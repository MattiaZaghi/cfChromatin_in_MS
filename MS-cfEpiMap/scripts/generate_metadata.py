#!/usr/bin/env python3
"""
Generate config/sample_metadata.tsv from SNAPIE fragment BED files.

Steps performed automatically:
  1. Read all SNAPIE fragment BEDs from FRAGS_DIR.
  2. Apply QC filter: keep only samples where qc_pass=TRUE in QC_SUMMARY.
  3. For each patient, keep the best plasma sample (highest fragment count).
     If no plasma passes QC, fall back to the best CSF sample.
  4. Infer sex from chrY: if chrY fragments / total fragments > 1 % → M else F.
  5. Write config/sample_metadata.tsv — clinical columns left as NA.

Run once from the MS-cfEpiMap/ directory:
    python scripts/generate_metadata.py
"""

import os
import re
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
FRAGS_DIR    = "/date/gcb/gcb_MZ/SNAPIE/frags"
BAM_DIR      = "/date/gcb/gcb_MZ/SNAPIE/align/dac"
BAM_SUFFIX   = ".dac_filtered.dedup.unique.sorted.bam"
PEAKS_DIR    = "/date/gcb/gcb_MZ/SNAPIE/peaks"
PEAKS_SUFFIX = "_peaks.narrowPeak"
QC_SUMMARY   = "/date/gcb/gcb_MZ/SNAPIE/qc/qc_summary.tsv"
OUTPUT       = "config/sample_metadata.tsv"

CHRY_THRESHOLD = 0.0013  # fraction of fragments on chrY → Male
                         # Natural gap in this cohort: females cluster ≤0.00109, males ≥0.00146


def parse_sample_id(sid):
    """Decode group / sample_type / protocol / volume from a SNAPIE sample ID."""
    if sid == "MS-P-Mix":
        return dict(patient_id="Mix", sample_type="plasma",
                    group="QC_Mix", protocol="standard", volume_ul="NA")

    protocol  = "standard"
    volume_ul = "NA"
    working   = sid

    m = re.search(r'-(\d+)$', working)
    if m:
        volume_ul = m.group(1)
        working   = working[:m.start()]

    if working.endswith("-pA"):
        protocol = "pA";  working = working[:-3]
    elif working.endswith("-1D"):
        protocol = "1D";  working = working[:-3]
    else:
        m = re.search(r'-(V\d+)$', working)
        if m:
            protocol = m.group(1);  working = working[:m.start()]

    tokens     = working.split("-")
    patient_id = None
    sample_type = "unknown"
    group_str   = ""

    for i, tok in enumerate(tokens):
        if tok in ("P", "C") and i >= 1:
            rest = "-".join(tokens[i + 1:])
            if rest.startswith("MS-") or rest.startswith("Ctrl"):
                patient_id  = "-".join(tokens[:i])
                sample_type = "plasma" if tok == "P" else "CSF"
                group_str   = rest
                break

    if patient_id is None:
        for i in range(len(tokens)):
            rest = "-".join(tokens[i:])
            if rest.startswith("MS-") or rest == "Ctrl":
                patient_id = "-".join(tokens[:i]) if i > 0 else tokens[0]
                group_str  = rest
                break
        if patient_id is None:
            patient_id = sid;  group_str = sid

    if   "Ctrl"       in group_str:                                        group = "Ctrl"
    elif "MS-New"     in group_str:                                        group = "NEW"
    elif "Rituximab"  in group_str and ("Progressive" in group_str or "Prog" in group_str):
                                                                            group = "MS-Rituximab-Progressive"
    elif "Rituximab"  in group_str and "Stable" in group_str:              group = "MS-Rituximab-Stable"
    else:                                                                   group = "UNKNOWN"

    return dict(patient_id=patient_id, sample_type=sample_type,
                group=group, protocol=protocol, volume_ul=volume_ul)


def load_qc(path):
    """Return dict: sample_id → {n_fragments, qc_pass}."""
    qc = {}
    if not os.path.exists(path):
        print(f"WARNING: QC summary not found at {path}; all samples treated as passing.")
        return qc
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            qc[row["sample_id"]] = {
                "n_fragments": int(row["n_fragments"]) if row["n_fragments"].isdigit() else 0,
                "qc_pass":     row["qc_pass"].upper() == "TRUE",
            }
    return qc


def infer_sex(bed_path):
    """Return 'M' if chrY fraction > CHRY_THRESHOLD, else 'F'."""
    try:
        total = int(subprocess.check_output(
            f"wc -l < {bed_path}", shell=True).decode().strip())
        if total == 0:
            return "NA"
        chry = int(subprocess.check_output(
            f"grep -c '^chrY' {bed_path} || true", shell=True).decode().strip())
        return "M" if (chry / total) >= CHRY_THRESHOLD else "F"
    except Exception:
        return "NA"


COLUMNS = [
    "sample_id", "patient_id", "group", "sample_type", "protocol",
    "fragment_bed", "bam", "peaks",
    "batch", "input_volume_ul",
    "age", "sex", "disease_duration_yrs",
    "rituximab_months", "months_since_rtx",
    "edss", "nfl", "cd19_count",
    "qc_include",
]


def main():
    # ── 1. Discover all SNAPIE fragment BEDs ─────────────────────────────────
    all_samples = sorted(
        f.replace(".bed", "")
        for f in os.listdir(FRAGS_DIR) if f.endswith(".bed")
    )
    print(f"Found {len(all_samples)} fragment BEDs in {FRAGS_DIR}")

    # ── 2. Load QC results ────────────────────────────────────────────────────
    qc = load_qc(QC_SUMMARY)

    # ── 3. Parse metadata + apply QC filter ───────────────────────────────────
    parsed = []
    for sid in all_samples:
        info = parse_sample_id(sid)
        if info["group"] == "QC_Mix":
            continue
        qc_info   = qc.get(sid, {"n_fragments": 0, "qc_pass": True})
        if not qc_info["qc_pass"]:
            continue
        parsed.append({**info, "sample_id": sid,
                       "n_fragments": qc_info["n_fragments"]})

    print(f"QC-passing samples: {len(parsed)}")

    # ── 4. Pick best sample per patient (plasma first, then CSF) ─────────────
    by_patient = defaultdict(list)
    for s in parsed:
        by_patient[s["patient_id"]].append(s)

    selected = []
    for pid, samples in sorted(by_patient.items()):
        plasma = [s for s in samples if s["sample_type"] == "plasma"]
        csf    = [s for s in samples if s["sample_type"] == "CSF"]
        pool   = plasma if plasma else csf
        if not pool:
            print(f"  WARNING: patient {pid} has no QC-passing samples; skipping")
            continue
        best   = max(pool, key=lambda s: s["n_fragments"])
        selected.append(best)

    print(f"Selected {len(selected)} samples (one per patient, best by fragment count)")

    # ── 5. Infer sex from chrY ────────────────────────────────────────────────
    print("Inferring sex from chrY fragment fraction ...")
    rows = []
    for s in selected:
        sid      = s["sample_id"]
        bed_path = f"{FRAGS_DIR}/{sid}.bed"
        sex      = infer_sex(bed_path)
        rows.append({
            "sample_id":            sid,
            "patient_id":           s["patient_id"],
            "group":                s["group"],
            "sample_type":          s["sample_type"],
            "protocol":             s["protocol"],
            "fragment_bed":         bed_path,
            "bam":                  f"{BAM_DIR}/{sid}{BAM_SUFFIX}",
            "peaks":                f"{PEAKS_DIR}/{sid}{PEAKS_SUFFIX}",
            "batch":                "NA",
            "input_volume_ul":      s["volume_ul"],
            "age":                  "NA",
            "sex":                  sex,
            "disease_duration_yrs": "NA",
            "rituximab_months":     "NA",
            "months_since_rtx":     "NA",
            "edss":                 "NA",
            "nfl":                  "NA",
            "cd19_count":           "NA",
            "qc_include":           "TRUE",
        })
        print(f"  {sid}  →  sex={sex}  frags={s['n_fragments']:,}")

    # ── 6. Write TSV ──────────────────────────────────────────────────────────
    Path(OUTPUT).parent.mkdir(exist_ok=True)
    with open(OUTPUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    from collections import Counter
    groups = Counter(r["group"]  for r in rows)
    sexes  = Counter(r["sex"]    for r in rows)
    print(f"\nWritten {len(rows)} samples → {OUTPUT}")
    print("Groups:", dict(groups))
    print("Sex   :", dict(sexes))


if __name__ == "__main__":
    main()
