rule end_motif_gc:
    conda: "envs/align.yaml"
    input:
        bam=config['outputFolder'] + "/align/namesorted/{sample}.n_sorted.bam",
        genome=config.get('genome_fasta', config.get('refdir', 'ref_files') + "/fa/" + config.get('genome', 'hg19') + ".fa")
    output:
        motif=config['outputFolder'] + "/motifs/{sample}/{sample}_4NMER_bp_motif.bed"
    params:
        read_method=config.get('read_method', 'PE'),
        nmer=config.get('nmer', 4),
        outdir=config['outputFolder'] + "/motifs/{sample}/"
    shell:
        """
        mkdir -p {params.outdir}

        if [ "{params.read_method}" != "PE" ] || [ ! -s {input.bam} ] || [ ! -f "{input.genome}" ]; then
            touch {output.motif}
            exit 0
        fi

        SAMPLE="{wildcards.sample}"
        NMER="{params.nmer}"
        GENOME="{input.genome}"
        BAM="{input.bam}"
        OUTDIR="{params.outdir}"

        BEDPE="$OUTDIR/${{SAMPLE}}.bedpe"
        BEDFILT="$OUTDIR/${{SAMPLE}}_filtered.bedpe"
        BED_GC="$OUTDIR/${{SAMPLE}}_frags_gc.bed"
        BPR1="$OUTDIR/${{SAMPLE}}_${{NMER}}_bp_r1.bed"
        BPR2="$OUTDIR/${{SAMPLE}}_${{NMER}}_bp_r2.bed"
        BPR1FA="$OUTDIR/${{SAMPLE}}_${{NMER}}NMER_bp_r1.fa.bed"
        BPR2FA="$OUTDIR/${{SAMPLE}}_${{NMER}}NMER_bp_r2.fa.bed"

        # Generate BEDPE
        bedtools bamtobed -bedpe -i "$BAM" | \
            awk 'OFS="\t" {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $6-$2}}' | \
            awk '$11 >= 0' > "$BEDPE"

        # GC content
        awk 'OFS="\t" {{print $1, $2, $6, $7, $11}}' "$BEDPE" | \
            sort -k1,1 -k2,2n | \
            bedtools nuc -fi "$GENOME" -bed - | \
            awk 'OFS="\t" {{print $1, $2, $3, $4, $5, $7}}' > "$BED_GC"

        # Filter by GC
        awk 'BEGIN {{FS=OFS="\t"}} FNR==NR {{arr[$4]=$6;next}} ($7 in arr) {{print $0, arr[$7]}}' \
            "$BED_GC" "$BEDPE" > "$BEDFILT"

        # N-mer extraction
        awk -v nmer="$NMER" 'OFS="\t" {{print $1, $2, $2+nmer, $7, $8, $9, $11, $12, $1, $2, $6}}' "$BEDFILT" > "$BPR1"
        awk -v nmer="$NMER" 'OFS="\t" {{print $4, $6-nmer, $6, $7, $8, $10, $11, $12, $1, $2, $6}}' "$BEDFILT" > "$BPR2"

        paste "$BPR1" <(bedtools getfasta -fi "$GENOME" -bed "$BPR1" -s -tab -fo /dev/stdout | awk '{{print $2}}') | \
            awk 'OFS="\t" {{print $1, $10, $11, $6, $7, $8, toupper($12)}}' > "$BPR1FA"
        paste "$BPR2" <(bedtools getfasta -fi "$GENOME" -bed "$BPR2" -s -tab -fo /dev/stdout | awk '{{print $2}}') | \
            awk 'OFS="\t" {{print toupper($12)}}' > "$BPR2FA"

        awk '{{getline line < f2; print $0 "\t" line}}' f2="$BPR2FA" "$BPR1FA" > {output.motif}
        """
