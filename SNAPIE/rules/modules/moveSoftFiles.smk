rule moveSoftFiles_module:
    conda: "envs/common.yaml"
    output:
        done=config['outputFolder'] + "/modules/moveSoftFiles.done"
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo "moveSoftFiles placeholder" > {output.done}
        """
