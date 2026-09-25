"""
Module 1 — Build Reference Regulatory Element (RRE) Universe
=============================================================
Run once:  snakemake --cores 4 --use-conda build_rre_universe

Outputs committed to reference/rre/ and NOT re-run by the main pipeline.
Accession numbers and download URLs must be filled in SOURCES.md and in
scripts/module1_build_rre.py before executing this target.
"""

# If universe.bed / constitutive_anchors.bed (and, when configured,
# universe_top.bed) already exist on disk — e.g. built externally by
# build_tissue_markers_sadeh.py + universe_and_anchors.py for this profile —
# don't even define build_rre_universe. Skipping at rule-definition time
# (rather than inside the shell command) means the rule never enters the
# DAG: no job stats line, no SLURM submission, and no risk of a crash or
# interrupt deleting these pre-built files via Snakemake's on-failure output
# cleanup.
def _rre_file_ready(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0

_universe_top_bed = config["rre"].get("universe_top_bed", "")
_build_rre_universe_needed = not (
    _rre_file_ready(config["rre"]["universe_bed"])
    and _rre_file_ready(config["normalization"]["constitutive_anchors_bed"])
    and (not _universe_top_bed or _rre_file_ready(_universe_top_bed))
)

if _build_rre_universe_needed:
    rule build_rre_universe:
        """
        Download reference H3K27ac peaks, merge across cell types, apply filters,
        and produce all RRE subset BEDs.
        """
        input:
            housekeeping_seed = config["normalization"]["housekeeping_seed"],
            chrom_sizes       = config["genome"]["chrom_sizes"],
        output:
            universe          = config["rre"]["universe_bed"],
            cns               = config["rre"]["cns_bed"],
            immune            = config["rre"]["immune_bed"],
            bcell             = config["rre"]["bcell_bed"],
            gwas_proximal     = config["rre"]["gwas_proximal_bed"],
            constitutive      = config["normalization"]["constitutive_anchors_bed"],
            sources_log       = "reference/rre/SOURCES.md",
            # Per-cell-type CNS BEDs — empty placeholders for cell types not yet downloaded
            cns_celltypes      = expand(
                "reference/rre/{ct}_rre.bed",
                ct=config["rre"]["cns_cell_types"],
            ),
            # Per-cell-type B cell BEDs — exclusive signature regions per B cell type
            bcell_celltypes    = expand(
                "reference/rre/{ct}_rre.bed",
                ct=config["rre"]["bcell_cell_types"],
            ),
            # ChIPseeker-style annotation TSVs (inspection only, not used downstream)
            ann_universe       = "reference/rre/ms_rre_universe_annotated.tsv",
            ann_cns            = "reference/rre/cns_rre_annotated.tsv",
            ann_immune         = "reference/rre/immune_rre_annotated.tsv",
            ann_bcell          = "reference/rre/bcell_rre_annotated.tsv",
            ann_gwas           = "reference/rre/gwas_proximal_rre_annotated.tsv",
            ann_cns_celltypes  = expand(
                "reference/rre/{ct}_rre_annotated.tsv",
                ct=config["rre"]["cns_cell_types"],
            ),
            ann_bcell_celltypes = expand(
                "reference/rre/{ct}_rre_annotated.tsv",
                ct=config["rre"]["bcell_cell_types"],
            ),
        params:
            merge_dist  = config["rre"]["merge_distance_bp"],
            min_size    = config["rre"]["min_size_bp"],
            max_size    = config["rre"]["max_size_bp"],
            min_anchor  = config["normalization"]["min_anchor_celltypes"],
            dac_regions = config["genome"]["dac_regions"],
            outdir      = "reference/rre",
            # workflow.basedir is the repo root (where Snakefile lives), unaffected
            # by the workdir directive — this gives an absolute path to the script.
            script      = workflow.basedir + "/scripts/module1_build_rre.py",
        threads: config["threads"]
        resources:
            mem_mb = 16000,
        conda:  "../../envs/python_analysis.yaml"
        shell:
            """
            python {params.script} \
                --housekeeping_seed {input.housekeeping_seed} \
                --chrom_sizes       {input.chrom_sizes} \
                --dac_regions       {params.dac_regions} \
                --outdir            {params.outdir} \
                --merge_dist        {params.merge_dist} \
                --min_size          {params.min_size} \
                --max_size          {params.max_size} \
                --min_anchor_types  {params.min_anchor} \
                --threads           {threads}
            """
