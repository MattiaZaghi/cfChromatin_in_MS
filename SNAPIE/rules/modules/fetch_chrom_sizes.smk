import os as _os

_CS_REF_DIR = config.get('ref_dir', 'ref_files')
_CS_GENOME  = config.get('genome', 'hg19')

rule fetch_chrom_sizes:
    conda: "envs/common.yaml"
    output:
        chromsizes=_os.path.join(_CS_REF_DIR, f"{_CS_GENOME}.chrom.sizes")
    params:
        genome=_CS_GENOME,
        refdir=_CS_REF_DIR
    shell:
        """
        mkdir -p {params.refdir}
        if command -v fetchChromSizes >/dev/null 2>&1; then
            fetchChromSizes {params.genome} > {output.chromsizes} || true
        else
            echo "fetchChromSizes not available; creating empty sizes file" > {output.chromsizes}
        fi
        """
