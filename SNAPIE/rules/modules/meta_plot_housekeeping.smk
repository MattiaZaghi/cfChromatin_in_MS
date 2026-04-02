rule meta_plot_housekeeping:
    """
    Generate meta-coverage plots around housekeeping gene TSSs for each sample.
    Input: deeptools bigWig (normalized coverage).
    Output: per-sample PDF in reports/meta_plots/.
    """
    conda: "envs/r_plots.yaml"
    input:
        bw=config['outputFolder'] + "/bigwig/deeptools/{sample}.bw",
        regions=config.get('housekeeping_bed',
                           'housekeeping_genes_comprehensive.bed')
    output:
        pdf=config['outputFolder'] + "/reports/meta_plots/{sample}_housekeeping_meta.pdf"
    params:
        sample="{sample}",
        script=config.get('pathMetaPlotScript',
                          'auxiliar_programs/meta_plot_housekeeping.R'),
        metaplot_r=config.get('pathMetaPlotR',
                              'auxiliar_programs/meta_plot_r_functions.R'),
        window=config.get('meta_plot_window', 10000),
        binsize=config.get('meta_plot_binsize', 25),
        color=config.get('meta_plot_color', 'steelblue'),
        outdir=config['outputFolder'] + "/reports/meta_plots/"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.script} \
            --bw      {input.bw} \
            --regions {input.regions} \
            --output  {output.pdf} \
            --sample  {params.sample} \
            --window  {params.window} \
            --binsize {params.binsize} \
            --color   {params.color} \
            || true
        touch {output.pdf}
        """
