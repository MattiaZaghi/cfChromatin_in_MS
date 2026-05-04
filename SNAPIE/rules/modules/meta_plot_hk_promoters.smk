rule meta_plot_hk_promoters:
    """
    Per-sample meta-coverage plot at housekeeping gene promoters.
    Uses pre-centered 2 kb windows from housekeeping_promoters.bed.
    Style: black-background heatmap (sorted by signal) + red mean-profile,
    matching the Sadeh cfChIP-seq pipeline aesthetic.
    Output: reports/qc/hk_meta/{sample}_hk_promoters_meta.pdf
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        bw     = config['outputFolder'] + "/bigwig/{sample}.bw",
        bed    = config.get('housekeeping_promoters_bed',
                            'ref_files/housekeeping_promoters.bed'),
        script = "auxiliar_programs/meta_plot_hk_promoters.R"
    output:
        pdf = config['outputFolder'] + "/reports/qc/hk_meta/{sample}_hk_promoters_meta.pdf"
    params:
        sample = "{sample}"
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        Rscript {input.script} \
            --bw      {input.bw} \
            --bed     {input.bed} \
            --output  {output.pdf} \
            --sample  '{params.sample}' \
            --binsize 25 \
            2>&1 || true
        touch {output.pdf}
        """
