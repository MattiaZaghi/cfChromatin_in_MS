rule unique_frags:
    conda: "envs/common.yaml"
    input:
        bed=config['outputFolder'] + "/frags/{sample}.bed"
    output:
        csv=config['outputFolder'] + "/frags/{sample}/{sample}_unique_frags.csv"
    params:
        dirname=config['outputFolder'] + "/frags/{sample}"
    threads: 1
    shell:
        """
        mkdir -p {params.dirname}
        if [ -f {input.bed} ]; then \
            count=$(wc -l {input.bed} | cut -f1 -d' '); \
        else \
            count=0; \
        fi; \
        printf "sample,count\n{wildcards.sample},%s\n" "$count" > {output.csv}
        """
