rule quality_report_lite_module:
    conda: "envs/common.yaml"
    output:
        done=config['outputFolder'] + "/modules/quality_report_lite.done"
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo "quality_report_lite placeholder" > {output.done}
        """
