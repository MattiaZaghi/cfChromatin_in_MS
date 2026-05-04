# ── MultiQC — aggregate FastQC reports across all samples ─────────────────────
#
# This rule runs once, after ALL per-sample FastQC jobs have finished.
# It searches the fastqc/ output tree recursively so it picks up every
# per-sample directory regardless of internal filename conventions.

rule multiqc_fastqc:
    """Aggregate per-sample FastQC reports into a single MultiQC HTML report."""
    conda: "envs/multiqc.yaml"
    input:
        # Wait for every sample's FastQC sentinel before aggregating
        done=expand(
            config['outputFolder'] + "/fastqc/{sample}/.done",
            sample=SAMPLES
        )
    output:
        report=config['outputFolder'] + "/reports/multiqc/multiqc_fastqc.html"
    params:
        search_dir=config['outputFolder'] + "/fastqc",
        outdir=config['outputFolder'] + "/reports/multiqc",
        report_name="multiqc_fastqc"
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {params.search_dir} \
            --outdir {params.outdir} \
            --filename {params.report_name} \
            --force
        """
