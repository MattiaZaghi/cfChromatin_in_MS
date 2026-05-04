"""
Module 9 — Pathway Enrichment (GSEApy pre-ranked).
Snakemake script: called via script: directive.

PURPOSE
Translate the list of differentially enriched regulatory regions from
Module 7c into biological pathways. Because DESeq2 output is region-
level (RREs or bins), we first map each region to its nearest gene, then
rank genes by a signed significance metric and run GSEA.

APPROACH — Pre-ranked GSEA
Rather than a simple over-representation test on a binary list of
significant genes, pre-ranked GSEA uses the full ranked list, which
captures the continuous gradient of enrichment/depletion and avoids
arbitrary significance cut-offs for gene set membership.

Rank metric: signed(-log10(padj)) × sign(log2FC)
  Highly significant up-regulated regions → large positive score.
  Highly significant down-regulated regions → large negative score.
  Non-significant regions → near zero.

Gene sets tested:
  MSigDB Hallmarks (H)   — broad biological states
  MSigDB C7             — immunological signatures
  KEGG                  — curated pathways
  GO Biological Process — Gene Ontology terms

Results are plotted as a dot-plot (top N pathways sorted by NES, dot size
= number of matched genes, colour = -log10 FDR). Pathways with FDR < 0.25
are prioritised; if none meet that threshold, the top N by significance
are shown regardless.

NOTE: If a 'gene_name' column is absent from the DESeq2 results (which
requires a prior annotation step to map regions to nearest RefSeq gene),
the region_id string is used as a proxy. Pre-ranked GSEA will then only
find enrichment when gene-set gene names happen to match region IDs, which
is unlikely. Ensure annotation is run before this module for meaningful results.
"""

import warnings
warnings.filterwarnings("ignore")

import io
import os
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Step 1: Load DESeq2 results for this contrast ────────────────────────────
res      = pd.read_csv(snakemake.input.deseq2_res, sep="\t", index_col=0)
contrast  = snakemake.params.contrast
gene_sets = snakemake.params.gene_sets
top_n     = int(snakemake.params.top_n)
refseq_bed = getattr(snakemake.input, "refseq_genes", None)
# Also check the standard cache location
_cache_path = Path("reference/genome/refseq_tss_hg19.bed")
if refseq_bed is None and _cache_path.exists():
    refseq_bed = str(_cache_path)

# ── Step 2: Annotate regions with nearest RefSeq gene ─────────────────────────
# region_id format: "chr1:123-456"  (produced by DESeq2 from the count matrix index)
# We parse coordinates, write a temp BED, use bedtools closest to find the nearest
# RefSeq gene TSS, and attach the gene symbol.  If no annotation file is provided
# or bedtools fails, we fall back to the region_id as gene proxy (GSEA will be
# empty but will not crash).

def _fetch_refseq_bed_ucsc(cache_path: Path) -> str:
    """Download hg19 RefSeq gene TSS BED from UCSC if not cached."""
    if cache_path.exists():
        return str(cache_path)
    import urllib.request
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
    )
    print(f"[Module 9] Downloading RefSeq hg19 gene table from UCSC ...")
    gz_tmp = str(cache_path) + ".gz"
    urllib.request.urlretrieve(url, gz_tmp)
    import gzip
    rows = []
    with gzip.open(gz_tmp, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            # refGene columns: bin, name, chrom, strand, txStart, txEnd, ..., name2
            if len(parts) < 13:
                continue
            chrom, strand, txStart, txEnd = parts[2], parts[3], int(parts[4]), int(parts[5])
            gene_name = parts[12]
            tss = txStart if strand == "+" else txEnd
            rows.append(f"{chrom}\t{tss}\t{tss+1}\t{gene_name}\t0\t{strand}")
    os.remove(gz_tmp)
    # deduplicate by gene name, keep first occurrence
    seen = set()
    with open(cache_path, "w") as fh:
        for row in rows:
            gn = row.split("\t")[3]
            if gn not in seen:
                seen.add(gn)
                fh.write(row + "\n")
    subprocess.run(f"bedtools sort -i {cache_path} > {cache_path}.sorted && mv {cache_path}.sorted {cache_path}",
                   shell=True, check=True)
    print(f"[Module 9] Saved {len(seen)} RefSeq genes to {cache_path}")
    return str(cache_path)


def _annotate_regions(res_df: pd.DataFrame, genes_bed: str) -> pd.DataFrame:
    """Use bedtools closest to add gene_name column to res_df."""
    # Parse region_id: supports "chr1:123-456" or "chr1_123_456"
    def parse_region(rid):
        rid = str(rid)
        if ":" in rid and "-" in rid:
            chrom, rest = rid.split(":", 1)
            start, end  = rest.split("-")
            return chrom, int(start), int(end)
        parts = rid.replace(":", "_").split("_")
        if len(parts) >= 3:
            try:
                return parts[0], int(parts[1]), int(parts[2])
            except ValueError:
                pass
        return None, None, None

    rows_bed = []
    for rid in res_df.index:
        chrom, start, end = parse_region(rid)
        if chrom:
            rows_bed.append(f"{chrom}\t{start}\t{end}\t{rid}")

    if not rows_bed:
        return res_df

    tmp_regions = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    tmp_regions.write("\n".join(rows_bed) + "\n")
    tmp_regions.close()

    sorted_regions = tmp_regions.name + ".sorted"
    subprocess.run(f"bedtools sort -i {tmp_regions.name} > {sorted_regions}",
                   shell=True, check=True)

    result = subprocess.run(
        f"bedtools closest -a {sorted_regions} -b {genes_bed} -D b -t first",
        shell=True, capture_output=True, text=True,
    )
    os.remove(tmp_regions.name)
    os.remove(sorted_regions)

    rid_to_gene: dict = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 8:
            rid_to_gene[parts[3]] = parts[7]  # col 8 = gene_name from genes_bed

    res_df = res_df.copy()
    res_df["gene_name"] = res_df.index.map(rid_to_gene)
    return res_df


# Determine gene annotation source
if refseq_bed and os.path.exists(str(refseq_bed)):
    genes_bed_path = str(refseq_bed)
else:
    cache = Path("reference/genome/refseq_tss_hg19.bed")
    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        genes_bed_path = _fetch_refseq_bed_ucsc(cache)
    except Exception as e:
        print(f"[Module 9] WARNING: Could not obtain RefSeq BED: {e}. GSEA will use region IDs.")
        genes_bed_path = None

if "gene_name" not in res.columns:
    if genes_bed_path:
        print(f"[Module 9] Annotating {len(res):,} regions with nearest gene ...")
        res = _annotate_regions(res, genes_bed_path)
        n_annotated = res["gene_name"].notna().sum()
        print(f"[Module 9] Annotated {n_annotated:,}/{len(res):,} regions")
    else:
        res["gene_name"] = None

if "gene_name" in res.columns and res["gene_name"].notna().any():
    fallback = pd.Series(res.index.astype(str).values, index=res.index)
    res["gene"] = res["gene_name"].fillna(fallback)
else:
    print("[Module 9] WARNING: No gene annotation — GSEA will not match any gene sets.")
    res["gene"] = res.index.astype(str)

# ── Step 3: Compute the ranking metric ───────────────────────────────────────
# Drop regions with missing padj or log2FC (these are features where DESeq2
# did not produce a valid estimate, typically due to all-zero counts or outliers).
res = res.dropna(subset=["padj", "log2FoldChange"])

# Rank metric: sign of fold change × magnitude of significance.
# clip(lower=1e-300) prevents log10(0) = -inf for perfectly significant features.
res["rank_metric"] = (
    -np.log10(res["padj"].clip(lower=1e-300)) *
    np.sign(res["log2FoldChange"])
)

# ── Step 4: Collapse to unique gene names ─────────────────────────────────────
# Multiple RREs can be assigned to the same nearest gene. We keep the maximum
# rank metric per gene, which represents the most significant regulatory signal
# linked to that gene. This is conservative — the gene is scored by its
# strongest-evidence regulatory element.
ranked = (
    res.groupby("gene")["rank_metric"]
    .max()
    .sort_values(ascending=False)
)
print(f"[Module 9] {contrast}: {len(ranked)} ranked genes/regions")

# ── Step 5: Run pre-ranked GSEA against each gene set library ─────────────────
# gseapy.prerank wraps the standard GSEA algorithm.
#   rnk:             the pre-sorted rank Series (index = gene name, value = metric)
#   gene_sets:       library name string (fetched from MSigDB by gseapy)
#   permutation_num: number of permutations for null distribution (1000 is standard)
#   seed:            fixed for reproducibility
# Errors are caught per library so that one failing library does not abort
# analysis of the others (e.g. network issues, empty gene set overlaps).
Path(snakemake.output.dotplot).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.output.table).parent.mkdir(parents=True, exist_ok=True)

def _normalise_gsea_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Harmonise gseapy res2d column names across package versions."""
    rename = {
        # gseapy >= 0.10
        "FDR q-val": "fdr", "NOM p-val": "pvalue", "Tag %": "tag_pct",
        "Gene %": "gene_pct", "Lead_genes": "lead_genes",
        # gseapy < 0.10
        "pval": "pvalue", "fdr": "fdr", "matched_size": "matched_size",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    # derive matched_size from Tag % if absent
    if "matched_size" not in df.columns and "tag_pct" in df.columns:
        df["matched_size"] = (df["tag_pct"].str.rstrip("%").astype(float)
                              * df.get("geneset_size", 100) / 100).round().astype(int)
    if "matched_size" not in df.columns:
        df["matched_size"] = 10  # fallback constant for dot sizing
    # ensure Term column exists
    if "Term" not in df.columns and df.index.name == "Term":
        df = df.reset_index()
    return df


all_results = []
for gs in gene_sets:
    try:
        pre_res = gp.prerank(
            rnk             = ranked,
            gene_sets       = gs,
            processes       = snakemake.threads,
            permutation_num = 1000,
            outdir          = None,
            seed            = 42,
            verbose         = False,
        )
        df = _normalise_gsea_cols(pre_res.res2d.copy())
        df["gene_set_library"] = gs
        print(f"[Module 9] {gs}: {len(df)} gene sets tested, "
              f"{(df['fdr'] < 0.25).sum()} with FDR<0.25")
        all_results.append(df)
    except Exception as e:
        print(f"[Module 9] WARNING: {gs} failed — {e}")

# ── Step 6: Aggregate results and write the combined table ───────────────────
if not all_results:
    # All gene set libraries failed (e.g. no gene name annotation, network error).
    # Write empty outputs so the Snakemake target file is still created.
    pd.DataFrame().to_csv(snakemake.output.table, sep="\t")
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, "No GSEA results", ha="center", va="center")
    fig.savefig(snakemake.output.dotplot)
    plt.close(fig)
else:
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(snakemake.output.table, sep="\t", index=False)

    # ── Step 7: Dotplot of top N pathways ────────────────────────────────────
    # Prioritise pathways with FDR < 0.25 (the conventional GSEA threshold).
    # Dot size = number of genes from the ranked list matched to the gene set
    #   (capped at 200 for visual clarity).
    # Dot colour = -log10(FDR) — brighter = more significant.
    # x-axis = NES (Normalised Enrichment Score): positive = up-regulated
    #   pathway, negative = down-regulated.
    combined["neg_log_fdr"] = -np.log10(combined["fdr"].clip(lower=1e-10))
    top = (
        combined[combined["fdr"] < 0.25]
        .nlargest(top_n, "neg_log_fdr")
        .sort_values("NES", ascending=True)   # sort so the plot reads top-to-bottom
    )

    # Fall back to the top N most significant if nothing meets FDR < 0.25
    if top.empty:
        top = combined.nlargest(top_n, "neg_log_fdr").sort_values("NES", ascending=True)

    fig, ax = plt.subplots(figsize=(7, max(4, len(top) * 0.32)))
    scatter = ax.scatter(
        top["NES"],
        range(len(top)),
        c  = top["neg_log_fdr"],
        s  = top["matched_size"].clip(upper=200) * 2,   # area ∝ matched gene set size
        cmap = "YlOrRd",
        vmin = 0, vmax = top["neg_log_fdr"].max(),
    )
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["Term"], fontsize=8)
    ax.axvline(0, ls="--", lw=0.8, color="grey")
    ax.set_xlabel("Normalised Enrichment Score (NES)")
    ax.set_title(f"GSEA — {contrast}  (top {len(top)} pathways)")
    plt.colorbar(scatter, ax=ax, label="-log10(FDR)")
    plt.tight_layout()
    fig.savefig(snakemake.output.dotplot, dpi=150)
    plt.close(fig)

print(f"[Module 9] {contrast} done.")
