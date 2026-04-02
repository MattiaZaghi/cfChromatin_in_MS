rule snp_smash_fingerprint_module:
    conda: "envs/common.yaml"
    output:
        done=config['outputFolder'] + "/modules/snp_smash_fingerprint.done"
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo "snp_smash_fingerprint placeholder" > {output.done}
        """
