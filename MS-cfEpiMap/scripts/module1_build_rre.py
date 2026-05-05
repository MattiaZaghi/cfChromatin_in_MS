#!/usr/bin/env python3
"""
Module 1 — Build Reference Regulatory Element (RRE) Universe.

Downloads H3K27ac ChIP-seq peak files from BLUEPRINT, ENCODE GRCh37, and
Roadmap Epigenomics, merges them across cell types, applies filters, and
produces all RRE subset BEDs.

BEFORE RUNNING: Fill in the accession numbers and download URLs in the
REFERENCE_DATA dict below. See reference/rre/SOURCES.md for the full list.

Usage (called by Snakemake module1_rre_universe.smk):
    python scripts/module1_build_rre.py \\
        --housekeeping_seed reference/rre/housekeeping_seed.bed \\
        --chrom_sizes /home/mattia/Genomes/hg19/hg19.chrom.sizes \\
        --dac_regions /home/mattia/Genomes/hg19/hg19_dac_exclusion.bed \\
        --outdir reference/rre \\
        --merge_dist 200 --min_size 100 --max_size 10000 \\
        --min_anchor_types 14 --threads 8

WHY BUILD A REFERENCE UNIVERSE?
Each individual cfChIP-seq sample typically yields fewer than 1,000 MACS2 peaks
due to the low input of cell-free DNA. This is too sparse for peak-centric
differential analysis with DESeq2 (which needs well-powered count distributions).
Instead, we define the feature universe from published H3K27ac ChIP-seq data of
the 15–16 relevant cell types (CNS, immune) using bulk tissue with high signal.
This gives ~300K–600K consistently defined regulatory elements that we then
COUNT our cfChIP-seq fragments against, giving robust integer matrices for DESeq2.
"""

import argparse
import os
import subprocess
import tempfile
import urllib.request
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# FILL IN ACCESSION NUMBERS / URLs BEFORE RUNNING
# Key: cell-type label used for annotation  Value: dict with url and tags
#
# All URLs must point to hg19-coordinate peak files (BED or narrowPeak, gzipped).
# BLUEPRINT and Roadmap Epigenomics are natively hg19.
# ENCODE data should be downloaded as GRCh37 (= hg19) versions.
# Neuro-epigenomics datasets originally in hg38 require prior liftOver with
# hg38ToHg19.over.chain before being placed here.
# ─────────────────────────────────────────────────────────────────────────────
REFERENCE_DATA = {
    # ── BLUEPRINT (hg19 native) ───────────────────────────────────────────
    "B_naive":        {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E032-H3K27ac.narrowPeak.gz",
                       "tags": ["bcell", "immune"]},
    "B_memory":       {"url": "TODO_BLUEPRINT_B_memory_H3K27ac_peaks.bed.gz",
                       "tags": ["bcell", "immune"]},
    "CD4_Th1":        {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E041-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},
    "CD4_Th17":       {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E042-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},
    "CD4_Treg":       {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E044-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},
    "CD8_T_naive":    {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E047-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},
    "CD8_T_memory":   {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E048-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},                   
    "Monocyte_CD14":  {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E029-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},
    "Monocyte_CD16":  {"url": "TODO_BLUEPRINT_CD16mono_H3K27ac_peaks.bed.gz",
                       "tags": ["immune"]},
    "NK":             {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E046-H3K27ac.narrowPeak.gz",
                       "tags": ["immune"]},
    "Neutrophil":     {"url": "TODO_BLUEPRINT_Neutrophil_H3K27ac_peaks.bed.gz",
                       "tags": ["immune"]},
    "Megakaryocyte":  {"url": "TODO_BLUEPRINT_Megakaryocyte_H3K27ac_peaks.bed.gz",
                       "tags": ["other"]},
    "PBMC;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E062-H3K27ac.narrowPeak.gz",
                        "tags": ["immune"]},
    "T-Cells;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E047-H3K27ac.narrowPeak.gz",
                                                  "tags": ["immune"]},
    "Treg-Cells;T-Cells;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E044-H3K27ac.narrowPeak.gz",
                                                             "tags": ["immune"]},
    "T-Helper-Cells;T-Cells;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E038-H3K27ac.narrowPeak.gz",
                                                                 "tags": ["immune"]},
    "T-Helper-Cells;T-Memory-Cells;T-Cells;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E037-H3K27ac.narrowPeak.gz",
                                                                                 "tags": ["immune"]},
    "T-Memory-Cells;T-Cells;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E048-H3K27ac.narrowPeak.gz",
                                                                  "tags": ["immune"]},
    "Monocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E029-H3K27ac.narrowPeak.gz",
                              "tags": ["immune"]},
    "B-Cells;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E032-H3K27ac.narrowPeak.gz",
                                                   "tags": ["immune"]},
    "NK;Lymphocytes;Leukocytes;Leukocytes": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E046-H3K27ac.narrowPeak.gz",
                                               "tags": ["immune"]},
    "Adipose": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E063-H3K27ac.narrowPeak.gz",
                "tags": ["other"]},
    "Atrium;Heart": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E104-H3K27ac.narrowPeak.gz",
                      "tags": ["other"]},
    "Ventricle;Heart": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E105-H3K27ac.narrowPeak.gz",
                         "tags": ["other"]},
    "Aorta;Heart": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E065-H3K27ac.narrowPeak.gz",
                     "tags": ["other"]},
    "GI Sm. Muscle": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E111-H3K27ac.narrowPeak.gz",
                       "tags": ["other"]},
    "Colon;Digestive": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E106-H3K27ac.narrowPeak.gz",
                         "tags": ["other"]},
    "Rectum;GI Mucosa;Digestive": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E102-H3K27ac.narrowPeak.gz",
                                    "tags": ["other"]},
    "Stomach;Digestive": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E094-H3K27ac.narrowPeak.gz",
                           "tags": ["other"]},
    "Placenta": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E091-H3K27ac.narrowPeak.gz",
                 "tags": ["other"]},
    "Ovary": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E097-H3K27ac.narrowPeak.gz",
              "tags": ["other"]},
    "Pancreas Islet;Pancreas": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E087-H3K27ac.narrowPeak.gz",
                                "tags": ["other"]},
    "Liver": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E066-H3K27ac.narrowPeak.gz",
              "tags": ["other"]},
    "Pancreas": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E098-H3K27ac.narrowPeak.gz",
                 "tags": ["other"]},
    "Lung": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E128-H3K27ac.narrowPeak.gz",
             "tags": ["other"]},
    "Epithelial": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E119-H3K27ac.narrowPeak.gz",
                    "tags": ["other"]},
    "Vasculary": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E122-H3K27ac.narrowPeak.gz",
                   "tags": ["other"]},
    "Skin": {"url": "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E127-H3K27ac.narrowPeak.gz",
             "tags": ["other"]},
    # ── Anita single-cell (CNS) ───────────────────────────────────────────────
    "Astrocyte":{"url": "/Users/gcblab/Anita/all_markers/AST.bed",
                       "tags": ["cns"]},
    "FIB":            {"url": "/Users/gcblab/Anita/all_markers/FIB.bed",
                       "tags": ["cns"]},
    "Microglia":         {"url": "/Users/gcblab/Anita/all_markers/MIGL-PVM.bed",
                       "tags": ["cns"]},
    "Oligodendrocyte":      {"url": "/Users/gcblab/Anita/all_markers/MOL.bed",
                       "tags": ["cns"]},
    "Neuron":      {"url": "/Users/gcblab/Anita/all_markers/NEU.bed",
                       "tags": ["cns"]},
    "OPC":      {"url": "/Users/gcblab/Anita/all_markers/OPC.bed",
                       "tags": ["cns"]},
    "VLMC":      {"url": "/Users/gcblab/Anita/all_markers/VLMC.bed",
                       "tags": ["cns"]}
}
# ─────────────────────────────────────────────────────────────────────────────


def run(cmd, check=True):
    """Run a shell command; raise RuntimeError on non-zero exit if check=True."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\n{result.stderr}")
    return result.stdout


def download_peaks(url: str, dest: str) -> bool:
    """
    Fetch a peak file to dest.
    Accepts:
      - http/https/ftp URLs  → downloaded with urllib
      - Local absolute paths → copied with shutil
    Returns False (and skips) if the URL is still a TODO placeholder.
    """
    if url.startswith("TODO"):
        print(f"  SKIP (TODO placeholder): {url}")
        return False
    if os.path.isfile(url):
        # In-house or pre-downloaded file: just copy it into raw_dir
        import shutil
        shutil.copy(url, dest)
        print(f"  Copied local file: {url}")
        return True
    print(f"  Downloading {url}")
    urllib.request.urlretrieve(url, dest)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--housekeeping_seed", required=True,
                    help="BED of housekeeping loci to seed the constitutive anchor set")
    ap.add_argument("--chrom_sizes",       required=True,
                    help="hg19 chromosome sizes (for boundary checks)")
    ap.add_argument("--dac_regions",       required=True,
                    help="ENCODE DAC exclusion list in hg19")
    ap.add_argument("--outdir",            required=True,
                    help="Output directory (reference/rre/)")
    ap.add_argument("--merge_dist",  type=int, default=200,
                    help="Max gap (bp) between peaks to merge into one element")
    ap.add_argument("--min_size",    type=int, default=100,
                    help="Minimum element size after merging (bp)")
    ap.add_argument("--max_size",    type=int, default=10000,
                    help="Maximum element size after merging (bp); removes super-large artefacts")
    ap.add_argument("--min_anchor_types", type=float, default=0.8,
                    help="Constitutive anchor threshold. If < 1: fraction of downloaded "
                         "cell types that must be active (e.g. 0.8 = 80%%). "
                         "If >= 1: absolute count (capped to the number of downloaded types).")
    ap.add_argument("--threads",     type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw_dir = outdir / "raw"   # store individual downloaded files here
    raw_dir.mkdir(exist_ok=True)

    n_total     = len(REFERENCE_DATA)
    downloaded  = []     # list of (cell_type, local_file_path) for downloaded files
    ct_tag_map  = {}     # cell_type → list of biological tags (cns, immune, bcell, …)

    # ── Step 1: Download reference peak files ─────────────────────────────────
    # Each cell type's H3K27ac peak file is fetched from BLUEPRINT/ENCODE/Roadmap.
    # Files that still have TODO placeholder URLs are skipped with a warning.
    print(f"\n[Module 1] Downloading {n_total} reference peak sets …")
    for ct, info in REFERENCE_DATA.items():
        dest = raw_dir / f"{ct}_peaks.bed.gz"
        ok   = download_peaks(info["url"], str(dest))
        if ok:
            downloaded.append((ct, str(dest)))
            ct_tag_map[ct] = info["tags"]

    # If no files were downloaded yet (all TODO), write placeholder outputs
    # so downstream Snakemake rules can still be developed and tested.
    if not downloaded:
        print(
            "[Module 1] WARNING: No reference peaks downloaded (all are TODO placeholders). "
            "Fill in REFERENCE_DATA URLs in scripts/module1_build_rre.py before running. "
            "Using housekeeping seed as constitutive anchors placeholder."
        )
        _write_placeholder_outputs(args, outdir)
        return

    # ── Step 2: Tag each peak with its cell-type and merge across all types ───
    # For each cell type, the peak BED is decompressed and a 4th column of
    # comma-separated tags (cell-type name + biological group tags) is appended.
    # All tagged peaks are concatenated, sorted, and merged with bedtools merge.
    # -d merge_dist merges peaks within merge_dist bp of each other — this
    # bridges nearby peaks from different cell types into a single element.
    # -c 4 -o collapse retains all unique tag strings from merged peaks.
    print("\n[Module 1] Merging peaks …")
    all_tagged = []
    for ct, bed_gz in downloaded:
        tags_str = ",".join([ct] + ct_tag_map[ct])   # e.g. "B_naive,bcell,immune"
        tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False)
        # zcat decompresses; awk reformats to 4-column BED with the tag in col 4
        run(f"zcat {bed_gz} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,\"{tags_str}\"}}' "
            f">> {tmp.name}")
        all_tagged.append(tmp.name)

    # Concatenate all cell-type BEDs, sort genome-wide, and merge overlapping peaks
    combined = outdir / "_combined.bed"
    run(f"cat {' '.join(all_tagged)} | bedtools sort -i - > {combined}")
    for t in all_tagged:
        os.remove(t)   # clean up individual temp files

    merged_raw = outdir / "_merged_raw.bed"
    run(
        f"bedtools merge -i {combined} -d {args.merge_dist} "
        f"-c 4 -o collapse > {merged_raw}"
        # -c 4 -o collapse: collapse all tag strings from merged peaks into one field
    )
    os.remove(combined)

    # ── Step 3: Size filter and DAC exclusion ────────────────────────────────
    # Elements shorter than min_size are typically noise or narrow artefacts.
    # Elements longer than max_size are usually repeat regions or assembly issues.
    # DAC (Difficult-to-map Accessible Chromatin) exclusion regions are removed
    # because they have mappability artefacts that inflate counts systematically.
    filtered = outdir / "_filtered.bed"
    run(
        f"awk 'BEGIN{{OFS=\"\\t\"}} "
        f"{{len=$3-$2; if (len>={args.min_size} && len<={args.max_size}) print}}' "
        f"{merged_raw} > {filtered}"
    )

    # Apply DAC removal if the exclusion file exists
    dac_cmd = (
        f"bedtools subtract -a {filtered} -b {args.dac_regions} -A"
        if os.path.exists(args.dac_regions)
        else f"cat {filtered}"   # if no DAC file, pass through unchanged
    )
    universe = outdir / "_universe_tmp.bed"
    run(f"{dac_cmd} > {universe}")
    os.remove(merged_raw)
    os.remove(filtered)

    # ── Step 4: Derive biological subsets ────────────────────────────────────
    # The tag column (col 4) records which cell type(s) each merged element
    # originated from. Grepping for a tag extracts that subset of elements.
    # These subsets are EXCLUSIVE: a region must carry the target group tag
    # but must NOT carry any of the excluded group tags.  Regions shared across
    # groups (e.g. active in both CNS and immune cell types) are dropped from
    # both subsets — they carry mixed signal and would confound deconvolution.
    # They are still present in the full ms_rre_universe.bed.
    def extract_exclusive_subset(include_tag, exclude_tags, dest):
        """
        Keep regions with include_tag that have none of the exclude_tags.
        Piped grep -v calls act as sequential exclusion filters.
        """
        cmd = f"grep '{include_tag}' {universe}"
        for ex_tag in exclude_tags:
            cmd += f" | grep -v '{ex_tag}'"
        cmd += f" > {dest}"
        run(cmd, check=False)   # grep exits 1 when no lines match — that is fine
        n = sum(1 for _ in open(dest))
        excl_str = ", ".join(exclude_tags) if exclude_tags else "none"
        print(f"  {include_tag} (excluding: {excl_str}): {n} regions → {dest}")

    # CNS-specific: active in ≥1 CNS type, not active in any immune type
    extract_exclusive_subset("cns",    ["immune"],  str(outdir / "cns_rre.bed"))

    # Immune-specific: active in ≥1 immune type, not active in any CNS type
    extract_exclusive_subset("immune", ["cns"],     str(outdir / "immune_rre.bed"))

    # B cell-specific: active in B naive/memory, not active in any CNS type
    # (bcell regions can overlap with other immune types — that is expected)
    extract_exclusive_subset("bcell",  ["cns"],     str(outdir / "bcell_rre.bed"))

    # ── Step 4b: Per-cell-type CNS subsets ───────────────────────────────────
    # For each CNS cell type defined in REFERENCE_DATA, produce a dedicated BED
    # of regions exclusive to that cell type (active in it, not in immune cells).
    # This supports per-cell-type DESeq2 contrasts and fine-grained deconvolution.
    # Cell types whose URLs are still TODO get an empty placeholder file so all
    # Snakemake output declarations are always satisfied.
    print("\n[Module 1] Building per-cell-type CNS subsets …")
    all_cns_cts      = [ct for ct, info in REFERENCE_DATA.items()
                        if "cns" in info.get("tags", [])]
    downloaded_names = {ct for ct, _ in downloaded}

    for ct in all_cns_cts:
        dest = outdir / f"{ct.lower()}_rre.bed"
        if ct in downloaded_names:
            # Exclude ALL other downloaded cell types — both other CNS types and
            # immune types.  A region shared between e.g. Oligodendrocyte and
            # Neuron (or Oligodendrocyte and any immune cell) would confound
            # deconvolution signal, so it is removed from both files.
            # Only regions found exclusively in this one cell type are retained.
            all_others = [other for other in downloaded_names if other != ct]
            extract_exclusive_subset(ct, all_others, str(dest))
        else:
            dest.write_text("")
            print(f"  {ct}: 0 regions (no data downloaded) → {dest}")

    # GWAS-proximal subset: RREs within 50 kb of MS GWAS SNPs (IMSGC)
    # bedtools window with -w 50000 extends each SNP by 50 kb on both sides
    # and -u returns each RRE at most once if it overlaps any window.
    gwas_out  = outdir / "gwas_proximal_rre.bed"
    gwas_snps = outdir / "ms_gwas_snps_hg19.bed"
    if gwas_snps.exists():
        run(
            f"bedtools window -a {universe} -b {gwas_snps} -w 50000 -u "
            f"> {gwas_out}"
        )
    else:
        print(f"  GWAS SNP file not found — writing empty gwas_proximal_rre.bed")
        gwas_out.write_text("")

    # Copy the full filtered universe to its final destination
    import shutil
    shutil.copy(str(universe), str(outdir / "ms_rre_universe.bed"))
    os.remove(universe)

    # ── Step 5: Build constitutive anchor set ────────────────────────────────
    # Constitutive anchors are RREs active in a high fraction of downloaded cell
    # types. Because they are ubiquitously active, their fragment count per sample
    # reflects only library size — making them ideal internal normalisation controls.
    #
    # The threshold is fraction-based so it adapts automatically to however many
    # cell types were actually downloaded (i.e. not TODO placeholders).
    # If min_anchor_types >= 1 it is treated as an absolute count capped to
    # len(downloaded); if < 1 it is treated as a fraction of len(downloaded).
    n_downloaded = len(downloaded)
    if args.min_anchor_types < 1.0:
        effective_min = max(1, round(args.min_anchor_types * n_downloaded))
    else:
        effective_min = min(int(args.min_anchor_types), n_downloaded)
    print(
        f"\n[Module 1] Building constitutive anchors "
        f"(≥{effective_min}/{n_downloaded} cell types, "
        f"threshold={args.min_anchor_types}) …"
    )
    anchors_out = outdir / "constitutive_anchors.bed"

    run(
        f"awk -v min={effective_min} "
        f"'BEGIN{{OFS=\"\\t\"}} "
        f"{{n=split($4,a,\",\"); ct=0; seen=\",\"; "
        f"for(i=1;i<=n;i++){{if(index(seen,\",\"a[i]\",\")==0){{ct++;seen=seen a[i] \",\"}}}}; "
        f"if(ct>=min) print $1,$2,$3,$4}}' "
        f"{outdir}/ms_rre_universe.bed > {anchors_out}"
    )
    n_anchors = sum(1 for _ in open(anchors_out))
    print(f"  {n_anchors} constitutive anchor regions (before housekeeping seed merge)")

    # Append housekeeping seed loci from SNAPIE meta-plot reference, then
    # re-sort and merge to produce a clean, non-overlapping anchor BED.
    run(f"cat {args.housekeeping_seed} >> {anchors_out}")
    run(f"bedtools sort -i {anchors_out} | bedtools merge -i - > {anchors_out}.tmp")
    os.replace(f"{anchors_out}.tmp", str(anchors_out))

    # ── Step 6: Write SOURCES.md documenting all accession numbers ───────────
    # A record of exactly which files were used is important for reproducibility
    # and for Methods sections in the manuscript.
    sources = outdir / "SOURCES.md"
    with open(sources, "w") as fh:
        fh.write("# RRE Universe — Data Sources\n\n")
        fh.write("| Cell type | Tags | URL |\n|---|---|---|\n")
        for ct, info in REFERENCE_DATA.items():
            fh.write(f"| {ct} | {','.join(info['tags'])} | {info['url']} |\n")
        fh.write(f"\nGenerated: {__import__('datetime').date.today()}\n")

    print("\n[Module 1] Done. Outputs in reference/rre/")


def _write_placeholder_outputs(args, outdir):
    """
    Create empty placeholder BEDs by copying the housekeeping seed.
    This lets downstream Snakemake rules be developed and tested before
    the real reference peak accessions have been filled in.
    """
    import shutil
    seed = Path(args.housekeeping_seed)
    for name in [
        "ms_rre_universe.bed", "cns_rre.bed", "immune_rre.bed",
        "bcell_rre.bed", "gwas_proximal_rre.bed", "constitutive_anchors.bed",
    ]:
        dest = outdir / name
        shutil.copy(str(seed), str(dest))
        print(f"  Placeholder: {dest} (copied from housekeeping seed)")
    (outdir / "SOURCES.md").write_text(
        "# SOURCES — fill in URLs in scripts/module1_build_rre.py\n"
    )


if __name__ == "__main__":
    main()
