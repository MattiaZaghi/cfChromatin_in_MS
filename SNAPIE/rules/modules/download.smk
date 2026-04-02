import os as _os

_REF_DIR = config.get('ref_dir', 'ref_files')
_GENOME  = config.get('genome', 'hg19')

rule download_snps:
    conda: "envs/common.yaml"
    output:
        snps=_os.path.join(_REF_DIR, f"snps_{_GENOME}.vcf")
    params:
        snp=config.get('snp', ''),
        refdir=_REF_DIR
    shell:
        """
        mkdir -p {params.refdir}
        if echo "{params.snp}" | grep -qE '^https?://'; then
            wget -O {output.snps} "{params.snp}" || : > {output.snps}
        elif [ -f "{params.snp}" ]; then
            cp "{params.snp}" {output.snps}
        else
            : > {output.snps}
        fi
        """

rule download_tss:
    conda: "envs/common.yaml"
    output:
        tss=_os.path.join(_REF_DIR, f"tss_promoter_peaks_{_GENOME}.bed")
    params:
        tss=config.get('tssPromoterPeaks', ''),
        refdir=_REF_DIR
    shell:
        """
        mkdir -p {params.refdir}
        if echo "{params.tss}" | grep -qE '^https?://'; then
            wget -O {output.tss} "{params.tss}" || : > {output.tss}
        elif [ -f "{params.tss}" ]; then
            cp "{params.tss}" {output.tss}
        else
            : > {output.tss}
        fi
        """

rule download_dac:
    conda: "envs/common.yaml"
    output:
        dac=_os.path.join(_REF_DIR, f"{_GENOME}.DAC.bed")
    params:
        dac=config.get('dacList', ''),
        refdir=_REF_DIR
    shell:
        """
        mkdir -p {params.refdir}
        if echo "{params.dac}" | grep -qE '^https?://'; then
            tmpgz={params.refdir}/$(basename {params.dac}).gz
            wget -O "$tmpgz" "{params.dac}" && gunzip -f "$tmpgz" || : > {output.dac}
        elif [ -f "{params.dac}" ]; then
            cp "{params.dac}" {output.dac}
        else
            : > {output.dac}
        fi
        """

rule download_gtf:
    conda: "envs/common.yaml"
    output:
        gtf=_os.path.join(_REF_DIR, f"{_GENOME}.GeneAnotation.gtf")
    params:
        gtf=config.get('geneAnnotation', ''),
        refdir=_REF_DIR
    shell:
        """
        mkdir -p {params.refdir}
        if [ -f "{params.gtf}" ]; then
            cp "{params.gtf}" {output.gtf}
        elif echo "{params.gtf}" | grep -qE '^https?://'; then
            tmpgz={params.refdir}/$(basename {params.gtf}).gz
            wget -O "$tmpgz" "{params.gtf}" && gunzip -f "$tmpgz" || : > {output.gtf}
        else
            : > {output.gtf}
        fi
        """

rule download_genome:
    conda: "envs/common.yaml"
    output:
        fa=_os.path.join(_REF_DIR, f"{_GENOME}.fa")
    params:
        fa=config.get('faGZFile', ''),
        refdir=_REF_DIR
    shell:
        """
        mkdir -p {params.refdir}
        if [ -f "{params.fa}" ]; then
            cp "{params.fa}" {output.fa}
        elif echo "{params.fa}" | grep -qE '^https?://'; then
            tmpgz={params.refdir}/$(basename {params.fa}).gz
            wget -O "$tmpgz" "{params.fa}" && gunzip -f "$tmpgz" || : > {output.fa}
        else
            : > {output.fa}
        fi
        """
