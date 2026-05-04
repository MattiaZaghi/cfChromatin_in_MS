"""
Module 9 — Pathway Enrichment (GSEApy)
========================================
Pre-ranked GSEA on DESeq2 results for each contrast.
Gene sets: MSigDB Hallmarks, C7 immunological, KEGG, GO BP.
"""

rule run_gsea:
    """
    Rank genes by signed -log10(padj) × sign(log2FC) from the full RRE
    DESeq2 results, run pre-ranked GSEA, generate dotplot.
    """
    input:
        deseq2_res = "results/differential/full/{contrast}_all.tsv",
    output:
        dotplot    = "results/gsea/{contrast}_dotplot.pdf",
        table      = "results/gsea/{contrast}_gsea_results.tsv",
    params:
        contrast = "{contrast}",
        gene_sets = [
            "MSigDB_Hallmark_2020",
            "MSigDB_C7_Immunologic_Signature",
            "KEGG_2021_Human",
            "GO_Biological_Process_2021",
        ],
        top_n   = 15,
    conda:  "../../envs/python_analysis.yaml"
    threads: 4
    script: "../../scripts/module9_gsea.py"
