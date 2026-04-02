rule fragle_ct_estimation:
    conda: "envs/common.yaml"
    input:
        bams=expand(
            config['outputFolder'] + "/align/fragle/{sample}.filtered.fragle.bam",
            sample=config.get('_samples_', [])
        ),
        bais=expand(
            config['outputFolder'] + "/align/fragle/{sample}.filtered.fragle.bam.bai",
            sample=config.get('_samples_', [])
        )
    output:
        fragle=config['outputFolder'] + "/reports/fragle/Fragle.txt"
    params:
        outdir=config['outputFolder'] + "/reports/fragle/",
        fragle_indir=config['outputFolder'] + "/align/fragle/"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}

        # Rename .filtered.fragle.bam to plain .bam so Fragle can pick them up
        WORKDIR=$(mktemp -d)
        for bam in {params.fragle_indir}*.filtered.fragle.bam; do
            [ -f "$bam" ] || continue
            base=$(basename "$bam" .filtered.fragle.bam)
            cp "$bam" "$WORKDIR/$base.bam"
            bai="$bam.bai"
            [ -f "$bai" ] && cp "$bai" "$WORKDIR/$base.bam.bai"
        done

        if command -v python >/dev/null 2>&1 && [ -f /usr/src/app/main.py ]; then
            cd /usr/src/app
            python /usr/src/app/main.py \
                --input "$WORKDIR" \
                --output {params.outdir} \
                --mode R \
                --cpu {threads} \
                --threads {threads}
            mv {params.outdir}Fragle.csv {output.fragle} 2>/dev/null || \
                mv {params.outdir}Fragle.txt {output.fragle} 2>/dev/null || \
                touch {output.fragle}
        else
            echo "Fragle tool not available; creating empty output" > {output.fragle}
        fi

        rm -rf "$WORKDIR"
        """
