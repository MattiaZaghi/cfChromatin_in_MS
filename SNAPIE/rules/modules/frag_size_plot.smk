rule frag_size_persample:
    """
    Per-sample fragment size distribution plot.
    Shows normalized frequency vs fragment size with nucleosomal markers
    at 147, 294, 441 bp and a red dotted median line.
    Output: reports/qc/frag_size/{sample}_fragment_size.pdf
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        frag_txt=config['outputFolder'] + "/frags/{sample}/{sample}.fragment_sizes.txt",
        script="auxiliar_programs/frag_size_plot_persample.R"
    output:
        pdf=config['outputFolder'] + "/reports/qc/frag_size/{sample}_fragment_size.pdf"
    params:
        sample="{sample}"
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        Rscript {input.script} \
            {input.frag_txt} \
            {output.pdf} \
            '{params.sample}' \
            2>&1 || true
        touch {output.pdf}
        """


rule frag_size_plot_groups:
    """
    Group-level fragment size distribution summary.
    Plots mean ± 1 SD per condition (Ctrl / NEW / MS-Rituximab-Stable /
    MS-Rituximab-Progressive) with nucleosomal markers at 147, 294, 441 bp.
    Output: reports/qc/fragment_size_by_group.pdf
    """
    conda: "envs/r_plots_dplyr.yaml"
    input:
        frag_sizes=expand(
            config['outputFolder'] + "/frags/{sample}/{sample}.fragment_sizes.txt",
            sample=config.get('_samples_', [])
        ),
        script="auxiliar_programs/frag_size_plot.R"
    output:
        pdf=config['outputFolder'] + "/reports/qc/fragment_size_by_group.pdf"
    params:
        frags_dir=config['outputFolder'] + "/frags",
        samplesheet=config.get('samplesheet', '')
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        Rscript {input.script} \
            {params.frags_dir} \
            {output.pdf} \
            {params.samplesheet} \
            2>&1 || true
        touch {output.pdf}
        """
