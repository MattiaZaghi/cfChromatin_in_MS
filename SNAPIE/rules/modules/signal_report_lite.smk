rule signal_report_lite_module:
    conda: "envs/common.yaml"
    output:
        done=config['outputFolder'] + "/modules/signal_report_lite.done"
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo "signal_report_lite placeholder" > {output.done}
        """
