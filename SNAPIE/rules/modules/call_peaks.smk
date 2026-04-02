
rule call_peaks:
    conda: "envs/peaks.yaml"
    input:
        bam=config['outputFolder'] + "/align/dedup/{sample}.dedup.unique.sorted.bam"
    output:
        narrowPeak=config['outputFolder'] + "/peaks/{sample}.narrowPeak",
        xls=config['outputFolder'] + "/peaks/{sample}_peaks.xls"
    params:
        genome_size="hs",
        dirname=config['outputFolder'] + "/peaks/"
    shell:
        """
        echo "Running MACS2 for {wildcards.sample}";
        mkdir -p {params.dirname};
        macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g {params.genome_size} -f BAMPE -t {input.bam} -n {wildcards.sample} --bdg --outdir {params.dirname} || true; 
        touch {output.narrowPeak} {output.xls}
        """
