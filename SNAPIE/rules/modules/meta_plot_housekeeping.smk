rule meta_plot_housekeeping:
    """
    Generate meta-coverage plots around TSS and Enhancer regions for each sample.
    Matches the cfChIP-seq Sadeh pipeline style:
      - Left panel:  TSS meta-plot (Meta-genes.bed with CpG/expression groups)
      - Right panel: Enhancer meta-plot (Meta-enhancers.bed with tissue-breadth groups)
      - Black heatmap background, red color, line plot above.
    Input: deeptools bigWig (normalized coverage).
    Output: per-sample PDF in reports/meta_plots/.
    """
    conda: "envs/r_plots.yaml"
    input:
        bw=config['outputFolder'] + "/bigwig/{sample}.bw"
    output:
        pdf=config['outputFolder'] + "/reports/meta_plots/{sample}_housekeeping_meta.pdf"
    params:
        sample="{sample}",
        script=config.get('pathMetaPlotScript',
                          'auxiliar_programs/meta_plot_housekeeping.R'),
        metaplot_r=config.get('pathMetaPlotR', ''),
        meta_genes=config.get('meta_genes_bed', ''),
        meta_enhancers=config.get('meta_enhancers_bed', ''),
        binsize=config.get('meta_plot_binsize', 25),
        outdir=config['outputFolder'] + "/reports/meta_plots/"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        ARGS="--bw {input.bw} --output {output.pdf} --sample {params.sample} --binsize {params.binsize}"
        [ -n "{params.metaplot_r}" ]     && ARGS="$ARGS --metaplot_r {params.metaplot_r}"
        [ -n "{params.meta_genes}" ]     && ARGS="$ARGS --meta_genes {params.meta_genes}"
        [ -n "{params.meta_enhancers}" ] && ARGS="$ARGS --meta_enhancers {params.meta_enhancers}"
        Rscript {params.script} $ARGS || true
        touch {output.pdf}
        """
