#!/usr/bin/env python3
"""Export an Excel sheet of the MS samples with disease and treatment.

Controls are excluded by design.

There is no clinical metadata file in this project: disease and treatment are
encoded in the sample ID itself, so everything here is parsed from
`config/samplesheet.csv`. The ID layout is

    <code>-<P|C>-<disease block>-<protocol>-<input ng>[-<mark>]

e.g. 04061-C-MS-Rituximab-Stable-V1-100, Nat035-P-MS-Natalizumab-V4-1000.
Anything the ID does not record (EDSS, NfL, age, sex, months since infusion,
CD19 count) is simply not available and is left out rather than guessed.

Usage:
    export_ms_sample_table.py <samplesheet.csv> <qc_summary.tsv> <out.xlsx>
"""

import re
import sys

import pandas as pd

# Disease block -> (disease, treatment, disease course). Ordered longest-first
# so that MS-Rituximab-Progressive is not swallowed by a shorter prefix.
BLOCKS = [
    ("MS-Rituximab-Progressive", ("Multiple sclerosis", "Rituximab",      "Progressive")),
    ("MS-Rituximab-Stable",      ("Multiple sclerosis", "Rituximab",      "Stable")),
    ("MS-Natalizumab",           ("Multiple sclerosis", "Natalizumab",    "Not annotated")),
    ("MS-New-PPMS",              ("Multiple sclerosis", "Treatment-naive", "Primary progressive")),
    ("MS-New-RR",                ("Multiple sclerosis", "Treatment-naive", "Relapsing-remitting")),
    ("MS-Unknown",               ("Multiple sclerosis", "Not annotated",  "Not annotated")),
]

SAMPLE_TYPE = {"P": "Plasma", "C": "CSF"}


def parse(sample_id):
    """Return a metadata dict, or None if the sample is not MS."""
    block = None
    for key, vals in BLOCKS:
        if key in sample_id:
            block = vals
            break
    if block is None:
        return None  # control, or a disease block we do not know about

    disease, treatment, course = block

    code = sample_id.split("-", 1)[0]

    m = re.match(r"^[A-Za-z0-9_.]+-([PC])-", sample_id)
    stype = SAMPLE_TYPE.get(m.group(1), "Unknown") if m else "Unknown"

    m = re.search(r"V[0-9]|(?<=[^A-Za-z])[0-9][A-Z]", sample_id)
    protocol = m.group(0) if m else "Unknown"

    m = re.search(r"-(\d+)(?:-H3K4me3)?$", sample_id)
    input_ng = int(m.group(1)) if m else None

    return {
        "sample_code": code,
        "sample_id": sample_id,
        "disease": disease,
        "treatment": treatment,
        "disease_course": course,
        "sample_type": stype,
        "protocol": protocol,
        "input_ng": input_ng,
    }


def main():
    sheet_f, qc_f, out_f = sys.argv[1], sys.argv[2], sys.argv[3]

    sheet = pd.read_csv(sheet_f)
    qc = pd.read_csv(qc_f, sep="\t")

    recs = [r for r in (parse(s) for s in sheet["sampleId"]) if r is not None]
    ms = pd.DataFrame(recs)

    ms = ms.merge(
        sheet[["sampleId", "enrichment_mark"]].rename(columns={"sampleId": "sample_id"}),
        on="sample_id", how="left",
    )
    ms = ms.merge(
        qc[["sample_id", "n_fragments", "enrichment_score", "qc_pass"]],
        on="sample_id", how="left",
    )

    ms["qc_pass"] = ms["qc_pass"].map(
        lambda v: "PASS" if str(v).strip().upper() == "TRUE"
        else ("FAIL" if str(v).strip().upper() == "FALSE" else "not run")
    )
    ms["enrichment_score"] = ms["enrichment_score"].round(2)

    # sample_type sits up front with the other headline fields: it is one of the
    # first things anyone filters on, so it should not be buried behind the
    # disease columns.
    cols = ["sample_code", "sample_id", "sample_type", "disease", "treatment",
            "disease_course", "protocol", "input_ng", "enrichment_mark",
            "n_fragments", "enrichment_score", "qc_pass"]
    ms = ms[cols].sort_values(["treatment", "disease_course", "sample_code"]) \
                 .reset_index(drop=True)

    v4 = ms[ms["protocol"] == "V4"].reset_index(drop=True)

    counts = (ms.groupby(["treatment", "disease_course", "protocol"])
                .size().reset_index(name="n_samples")
                .sort_values(["treatment", "disease_course", "protocol"]))

    with pd.ExcelWriter(out_f, engine="openpyxl") as xw:
        ms.to_excel(xw, sheet_name="MS samples", index=False)
        v4.to_excel(xw, sheet_name="MS samples V4 only", index=False)
        counts.to_excel(xw, sheet_name="Counts", index=False)

        for name, df in (("MS samples", ms), ("MS samples V4 only", v4),
                         ("Counts", counts)):
            ws = xw.sheets[name]
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
            for i, col in enumerate(df.columns, start=1):
                width = max(len(str(col)), *(len(str(v)) for v in df[col])) + 2
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(width, 34)

    print(f"wrote {out_f}")
    print(f"  MS samples: {len(ms)}  (V4: {len(v4)})")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
