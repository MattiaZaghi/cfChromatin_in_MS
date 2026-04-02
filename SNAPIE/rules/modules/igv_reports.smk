rule igv_consolidate_report:
    conda: "envs/common.yaml"
    input:
        bigwigs=expand(
            config['outputFolder'] + "/bigwig/{sample}.bw",
            sample=config.get('_samples_', [])
        )
    params:
        out_dir=config.get('outputFolder', '') + "/reports/multiqc",
        header=config.get('multiqc_housekeeping_header', ''),
        bigwig_dir=config.get('outputFolder', '') + "/bigwig"
    output:
        html=config['outputFolder'] + "/reports/multiqc/igv_housekeeping_genes_mqc.html"
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -f "{params.header}" ]; then
            cp "{params.header}" {output.html}
        else
            echo "<html><body><h1>IGV Reports</h1>" > {output.html}
        fi
        for file in {params.out_dir}/*_igv_housekeeping_genes_report.html; do
            [ -e "$file" ] || continue
            base=$(basename "$file")
            link_text=$(basename "$file" "_igv_housekeeping_genes_report.html")
            echo "<a href='$base' target='_blank' class='btn btn-primary'>$link_text</a>" >> {output.html}
        done
        echo "</body></html>" >> {output.html}
        """
