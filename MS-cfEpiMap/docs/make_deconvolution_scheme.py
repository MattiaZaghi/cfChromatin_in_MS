"""
Presentation scheme for Module 4 — within-sample cell-type deconvolution.
Generates a vector (PDF) + raster (PNG) flow diagram.

Run:  /home/mattia/miniconda3_n/envs/common/bin/python make_deconvolution_scheme.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

# ── palette ───────────────────────────────────────────────────────────────────
C_INPUT = "#4878CF"   # inputs (blue)
C_PREP  = "#6ACC65"   # region prep (green)
C_MODEL = "#D65F5F"   # statistical model (red)
C_CALL  = "#B47CC7"   # calls / effect size (purple)
C_OUT   = "#E8A33D"   # outputs (amber)
C_QC    = "#7F7F7F"   # QC (grey)
C_EDGE  = "#2B2B2B"

fig, ax = plt.subplots(figsize=(15.5, 9.6))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, text, color, fs=10.5, bold_first=True, text_color="white"):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.6",
        linewidth=1.4, edgecolor=C_EDGE, facecolor=color, alpha=0.95, zorder=2))
    if bold_first and "\n" in text:
        head, body = text.split("\n", 1)
        ax.text(x + w/2, y + h*0.74, head, ha="center", va="center",
                fontsize=fs+0.5, fontweight="bold", color=text_color, zorder=3)
        ax.text(x + w/2, y + h*0.34, body, ha="center", va="center",
                fontsize=fs-1.5, color=text_color, zorder=3)
    else:
        ax.text(x + w/2, y + h/2, text, ha="center", va="center",
                fontsize=fs, fontweight="bold" if bold_first else "normal",
                color=text_color, zorder=3)

def arrow(x1, y1, x2, y2, style="-|>", color=C_EDGE, lw=1.8, ls="-", rad=0.0):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle=style, mutation_scale=16,
        linewidth=lw, color=color, linestyle=ls,
        connectionstyle=f"arc3,rad={rad}", zorder=1))

# ── title ───────────────────────────────────────────────────────────────────
ax.text(50, 97.5, "Within-sample cell-type deconvolution (group-agnostic)",
        ha="center", va="center", fontsize=17, fontweight="bold")
ax.text(50, 93.6,
        "Per sample × per cell type: is H3K27ac cfDNA enriched at this cell type's "
        "signature, above its own local genomic background?",
        ha="center", va="center", fontsize=10.5, color="#444444", style="italic")

# ── ROW 1 — inputs ────────────────────────────────────────────────────────────
box( 2, 80, 22, 9,
     "Fragment BED  (1 / sample)\nmidpoint = (start+end)/2", C_INPUT, fs=10)
box(28, 80, 30, 9,
     "RRE universe  (334k regions)\ncell-type–exclusive H3K27ac signatures", C_INPUT, fs=10)
box(62, 80, 36, 9,
     "chromHMM off-target BED (220k)\nnever-active windows  →  null territory", C_INPUT, fs=10)

# ── ROW 2 — region preparation ──────────────────────────────────────────────
box(28, 65, 30, 9,
     "Signature regions per cell type\nexclusive to one type in the panel", C_PREP, fs=10)
box(62, 65, 36, 9,
     "Background tiles\nsubtract signatures → 5 kb windows", C_PREP, fs=10)
arrow(43, 80, 43, 74)            # RRE -> signatures
arrow(80, 80, 80, 74)            # off-target -> bg tiles

# ── ROW 3 — counting + GC ────────────────────────────────────────────────────
box(28, 50, 70, 9,
     "Count fragment midpoints per region   (bedtools intersect, parallel over samples)\n"
     "compute GC content per region (bedtools nuc)", C_PREP, fs=10)
arrow(13, 80, 30, 59, rad=-0.15)         # frags -> counting
arrow(43, 65, 43, 59)                    # signatures -> counting
arrow(80, 65, 80, 59)                    # bg tiles -> counting

# ── ROW 4 — background model + matching (two parallel green->red feeders) ─────
box( 2, 34, 45, 10,
     "Sadeh local background rate\n"
     "trimmed-mean density,  genome → chr → 2 Mb hierarchy", C_MODEL, fs=10)
box(53, 34, 45, 10,
     "GC × length matched background\n"
     "10 tiles / signature region drawn from matching bin", C_MODEL, fs=10)
arrow(50, 50, 24, 44, rad=0.12)          # counting -> sadeh
arrow(66, 50, 76, 44, rad=-0.10)         # counting -> matched bg
# per-sample dispersion note
ax.text(24.5, 31.2, "+ per-sample NB dispersion α  (from background variance)",
        ha="center", va="center", fontsize=8.6, color=C_MODEL, style="italic")

# ── ROW 5 — NB GLM (the test) ─────────────────────────────────────────────────
box(20, 16, 60, 11,
     "Negative-binomial GLM   (per sample × cell type)\n"
     r"counts ~ 1 + is_signature      offset = log(local expected count)",
     C_MODEL, fs=11)
ax.text(50, 18.4,
        r"$\beta_{1}$ = log-enrichment of signature over its local background  →  one-sided upper-tail $p$",
        ha="center", va="center", fontsize=9.2, color="white", zorder=3)
arrow(24, 34, 38, 27, rad=0.10)          # sadeh -> glm (offset)
arrow(76, 34, 62, 27, rad=-0.10)         # matched bg -> glm (comparison)

# ── ROW 6 — calls + effect size + contribution ───────────────────────────────
box( 2, 2, 30, 9,
     "BH-FDR within sample\nenriched = q < 0.05", C_CALL, fs=10)
box(35, 2, 30, 9,
     "log2 fold-enrichment\npooled obs / expected counts", C_CALL, fs=10)
box(68, 2, 30, 9,
     "Fractional contribution\nexcess share across cell types", C_CALL, fs=10)
arrow(40, 16, 17, 11, rad=0.12)
arrow(50, 16, 50, 11)
arrow(60, 16, 83, 11, rad=-0.12)

# ── QC sidebar ────────────────────────────────────────────────────────────────
ax.add_patch(FancyBboxPatch(
    (2, 50.0), 24, 8.6, boxstyle="round,pad=0.4,rounding_size=1.6",
    linewidth=1.3, edgecolor=C_EDGE, facecolor=C_QC, alpha=0.93, zorder=2))
ax.text(14, 56.3, "QC / null checks", ha="center", va="center",
        fontsize=10, fontweight="bold", color="white", zorder=3)
ax.text(14, 52.6, "permuted-signature null\n+ anti-signal (expect log2FE ≈ 0)",
        ha="center", va="center", fontsize=8.2, color="white", zorder=3)
arrow(26, 54, 28, 54.5, color=C_QC, lw=1.4, ls=(0, (4, 2)))

# ── legend ────────────────────────────────────────────────────────────────────
legend_items = [
    ("Input data", C_INPUT), ("Region preparation", C_PREP),
    ("Statistical model", C_MODEL), ("Per-sample call / score", C_CALL),
    ("QC", C_QC),
]
handles = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c,
                  markeredgecolor=C_EDGE, markersize=13, label=l)
           for l, c in legend_items]
ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.02),
          ncol=5, frameon=False, fontsize=9.5, handletextpad=0.4,
          columnspacing=1.4)

plt.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"within_sample_deconvolution_scheme.{ext}",
                dpi=200, bbox_inches="tight")
print("written: within_sample_deconvolution_scheme.pdf / .png")
