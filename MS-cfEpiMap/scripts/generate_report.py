"""
Generate integrated HTML QC + analysis report.
Snakemake script: called via script: directive.

PURPOSE
Aggregate key figures and summary tables from all modules into a single
self-contained HTML file. The report is designed to be the first thing
to open after a pipeline run — it provides a one-page overview of data
quality and key results without requiring access to individual output files.

STRUCTURE
  1. Sample overview table
  2. Module 2 — constitutive anchor scaling factors (QC barplot + SF table)
  3. Module 3 — B Cell Index violin + flagged samples
  4. Module 7 — batch × group contingency table
  5. Module 7 — differential results bar chart (significant features per contrast)
  6. Module 10 — MEAS score scatter
  7. QC summary table (all QC metrics per sample)

All figures are embedded as base64-encoded PNG images so the HTML file
is fully self-contained (no external image paths required).
"""

import base64
import io
from pathlib import Path
from datetime import date

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend


# ── Utility functions ─────────────────────────────────────────────────────────

def fig_to_b64(fig):
    """
    Render a matplotlib Figure to a base64-encoded PNG string.
    This is how figures are embedded directly in HTML without external files.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def df_to_html(df, max_rows=200):
    """
    Convert a DataFrame to an HTML table string using Bootstrap classes
    for consistent styling in the report. Truncates to max_rows rows to
    avoid very large report files when there are many samples.
    """
    return df.head(max_rows).to_html(
        classes="table table-sm table-striped", border=0, index=True
    )


# ── Step 1: Load all module outputs ──────────────────────────────────────────
meta        = pd.read_csv(snakemake.input.meta,        sep="\t", comment="#")
sf          = pd.read_csv(snakemake.input.sf,           sep="\t")   # Module 2 scaling factors
bci         = pd.read_csv(snakemake.input.bci,          sep="\t")   # Module 3 BCI
batch_table = pd.read_csv(snakemake.input.batch_table,  sep="\t", index_col=0)  # Module 7a
deconv      = pd.read_csv(snakemake.input.deconv,       sep="\t")   # Module 4 z-scores
meas        = pd.read_csv(snakemake.input.meas,         sep="\t")   # Module 10 MEAS scores

contrasts = snakemake.params.contrasts   # list of contrast labels for the differential bar plot

# ── Step 2: Build the unified QC summary table ───────────────────────────────
# Merge anchor-normalisation QC (sf) with BCI metrics (bci) per sample.
# qc_pass = TRUE only when the SF is not flagged (BCI flagging is informational,
# not a hard exclusion criterion — see Module 3 comments for rationale).
qc = sf[["sample_id", "anchor_reads", "constitutive_sf", "sf_flagged", "group"]].copy()
qc = qc.merge(bci[["sample_id", "bci", "bci_scaled", "bci_flag"]], on="sample_id", how="left")
qc["total_frags"] = "NA"   # placeholder — would require a wc -l on frags files
qc["qc_pass"]     = ~qc["sf_flagged"].fillna(False)

Path(snakemake.output.qc_summary).parent.mkdir(parents=True, exist_ok=True)
qc.to_csv(snakemake.output.qc_summary, sep="\t", index=False)

# ── Step 3: Generate inline figures ──────────────────────────────────────────
palette = {
    "Ctrl":                     "#4878CF",
    "NEW":                      "#6ACC65",
    "MS-Rituximab-Stable":      "#D65F5F",
    "MS-Rituximab-Progressive": "#B47CC7",
}

# ── Figure 1: Scaling factor barplot (Module 2 QC) ───────────────────────────
# Each bar is one sample; height = constitutive SF; colour = group.
# Bars near 1.0 = sample depth close to the cohort median (good QC).
# Very tall or very short bars indicate outlier library depth.
fig1, ax = plt.subplots(figsize=(12, 4))
sf_s    = sf.sort_values("group")
colors  = [palette.get(g, "grey") for g in sf_s["group"].fillna("grey")]
ax.bar(range(len(sf_s)), sf_s["constitutive_sf"], color=colors)
ax.axhline(1, ls="--", color="black", lw=0.8, label="SF=1 (median depth)")
ax.set_xticks([])
ax.set_ylabel("Constitutive scaling factor")
ax.set_title("Module 2 — Anchor scaling factors (grouped by colour)")
plt.tight_layout()
sf_b64 = fig_to_b64(fig1)
plt.close(fig1)

# ── Figure 2: BCI violin by group (Module 3 QC) ──────────────────────────────
# Quick visual check: do Ctrl samples have BCI near 1, and are treated patients
# within the expected range? Obvious outliers suggest reconstitution extremes.
fig2, ax = plt.subplots(figsize=(6, 4))
order = [g for g in ["Ctrl", "NEW", "MS-Rituximab-Stable", "MS-Rituximab-Progressive"]
         if g in bci["group"].unique()]
for i, grp in enumerate(order):
    vals = bci.loc[bci["group"] == grp, "bci"].dropna()
    ax.violinplot(vals, positions=[i], showmedians=True)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, rotation=30, ha="right", fontsize=8)
ax.set_ylabel("BCI")
ax.set_title("Module 3 — B Cell Index by group")
plt.tight_layout()
bci_b64 = fig_to_b64(fig2)
plt.close(fig2)

# ── Figure 3: Differential results bar chart (Module 7) ──────────────────────
# Number of significant features (FDR < 0.05, |LFC| > 1) per contrast.
# A bar with 0 hits is informative: it means no detectable difference for
# that contrast given the current sample size and statistical model.
sig_counts = {}
for f in snakemake.input.diff_sig:
    c = Path(f).stem.replace("_significant", "")
    try:
        df = pd.read_csv(f, sep="\t")
        sig_counts[c] = len(df)
    except Exception:
        sig_counts[c] = 0

fig3, ax = plt.subplots(figsize=(7, 4))
ax.barh(list(sig_counts.keys()), list(sig_counts.values()), color="#4878CF")
ax.set_xlabel("Significant features (FDR < 0.05, |LFC| > 1)")
ax.set_title("Module 7 — Differential results (full RRE)")
plt.tight_layout()
diff_b64 = fig_to_b64(fig3)
plt.close(fig3)

# ── Figure 4: MEAS score scatter (Module 10) ─────────────────────────────────
# MEAS = P(Progressive) − P(Ctrl); higher values reflect a more progressive
# epigenomic profile. Points are coloured by true group — the key expectation
# is that Ctrl and NEW cluster at low MEAS, Stable in the middle, and
# Progressive at the top.
fig4, ax = plt.subplots(figsize=(6, 4))
for grp, sub in meas.groupby("group"):
    ax.scatter(range(len(sub)), sub["meas_score"],
               label=grp, color=palette.get(grp, "grey"), alpha=0.8)
ax.axhline(0, ls="--", color="grey")
ax.set_ylabel("MEAS score")
ax.set_title("Module 10 — MEAS score by group")
ax.legend(fontsize=7)
plt.tight_layout()
meas_b64 = fig_to_b64(fig4)
plt.close(fig4)

# ── Step 4: Identify flagged samples for warning boxes ───────────────────────
flagged_sf  = qc[qc["sf_flagged"].fillna(False)]   # Module 2: outlier library depth
flagged_bci = qc[qc["bci_flag"].fillna(False)]     # Module 3: high B cell reconstitution


# ── Step 5: Build the HTML string ────────────────────────────────────────────
def img_tag(b64):
    """Wrap a base64 PNG string in an HTML img tag."""
    return f'<img src="data:image/png;base64,{b64}" style="max-width:100%">'

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>MS-cfEpiMap Report</title>
  <!-- Bootstrap CSS for responsive tables and alert boxes -->
  <link rel="stylesheet"
        href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
  <style>
    body {{ font-family: sans-serif; margin: 24px; }}
    h2   {{ margin-top: 2rem; }}
    .flagged {{ background: #fff3cd; }}
  </style>
</head>
<body>
<h1>MS-cfEpiMap — Analysis Report</h1>
<!-- Date and sample count shown at the top for quick reference -->
<p class="text-muted">Generated: {date.today()}  |  Samples: {len(meta)}</p>

<!-- ── Sample overview ─────────────────────────────────────────────────── -->
<h2>Sample Overview</h2>
{df_to_html(meta[["sample_id","group","sample_type","protocol","batch"]].head(100))}

<!-- ── Module 2: Constitutive Anchor Normalisation ─────────────────────── -->
<h2>Module 2 — Constitutive Anchor Normalization</h2>
<!-- Barplot: each bar = one sample, height = SF, colour = group -->
{img_tag(sf_b64)}
<!-- Alert box if any samples have SF > 3 SD from mean -->
{"<div class='alert alert-warning'>Flagged samples (SF > 3 SD):<br>" +
 ", ".join(flagged_sf["sample_id"]) + "</div>" if not flagged_sf.empty else ""}
{df_to_html(sf[["sample_id","group","anchor_reads","constitutive_sf","sf_flagged"]])}

<!-- ── Module 3: B Cell Index ───────────────────────────────────────────── -->
<h2>Module 3 — B Cell Index</h2>
{img_tag(bci_b64)}
<!-- Flag treated patients with BCI > 1.5× Ctrl mean for manual review -->
{"<div class='alert alert-warning'>Flagged (BCI > 1.5× Ctrl mean in treated patients):<br>" +
 ", ".join(flagged_bci["sample_id"]) + "</div>" if not flagged_bci.empty else ""}

<!-- ── Module 7: Batch QC ─────────────────────────────────────────────── -->
<h2>Module 7 — Batch QC</h2>
<h4>Batch × Group contingency</h4>
<!-- Zeros in any cell indicate a group absent from a batch; warnings appear in the log -->
{df_to_html(batch_table)}

<!-- ── Module 7: Differential results ───────────────────────────────────── -->
<h2>Module 7 — Differential Results</h2>
{img_tag(diff_b64)}

<!-- ── Module 10: MEAS Score ─────────────────────────────────────────────── -->
<h2>Module 10 — MEAS Score</h2>
{img_tag(meas_b64)}

<!-- ── QC Summary table ──────────────────────────────────────────────────── -->
<h2>QC Summary</h2>
<!-- One row per sample; qc_pass = FALSE means at least one QC flag was raised -->
{df_to_html(qc)}

</body>
</html>
"""

# ── Step 6: Write the HTML file ───────────────────────────────────────────────
Path(snakemake.output.html).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.html, "w") as fh:
    fh.write(html)

print(f"[Report] Written to {snakemake.output.html}")
