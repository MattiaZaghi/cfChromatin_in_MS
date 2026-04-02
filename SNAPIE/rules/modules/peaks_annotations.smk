rule peaks_annotations:
    conda: "envs/peaks.yaml"
    input:
        script=config.get('pathRGenomicAnnotation', 'auxiliar_programs/genomic_annotation.R'),
        peaks=expand(
            config['outputFolder'] + "/peaks/{sample}.narrowPeak",
            sample=config.get('_samples_', [])
        )
    params:
        peaks_dir=config.get('outputFolder', '') + "/peaks",
        out_dir=config.get('outputFolder', '') + "/reports/multiqc"
    output:
        done=config['outputFolder'] + "/reports/multiqc/peaks_annotations.done"
    threads: 1
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -f "{input.script}" ]; then
            Rscript {input.script} -i {params.peaks_dir} -o {params.out_dir} --force_barplot TRUE || true
        fi
        touch {output.done}
        """
