include: 'Snakefile_prep.smk'

rule all_preprocess:
    input:
        cellranger=['{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'.format(sample=sample,modality=modality,barcode=barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        bigwig_all=['{sample}/{modality}_{barcode}/bigwig/all_reads.bw'.format(sample=sample,modality=modality,barcode=barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        macs_broad=['{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample,modality=modality,barcode=barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        peaks_overlap=['{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'.format(sample=sample,modality=modality,barcode=barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        barcodes_sum=['{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'.format(sample=sample,modality=modality,barcode=barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        cell_pick=['{sample}/{modality}_{barcode}/cell_picking/metadata.csv'.format(sample=sample,modality=modality,barcode=barcodes_dict[sample][modality]) for sample in samples_list for modality in barcodes_dict[sample].keys()],
        bam_RNA=['{sample}/RNA_AAAAGGGG_H3K27ac/cellranger/outs/possorted_genome_bam.bam'.format(sample=sample) for sample in samples_list],
        #agg=['{merge}/{modality_n}_{barcode_n}/cellranger/outs/fragments.tsv.gz'.format(merge=merge,modality_n=modality,barcode_n=barcodes_gen[merge][modality]) for merge in params_list for modality in barcodes_gen[merge].keys()],
        #bw_merge_ATAC =['{merge}/ATAC_TATAGCCT/bigwig/all_reads.bw'.format(merge=merge) for merge in params_list],
        #bw_merge_H3K27ac =['{merge}/H3K27ac_CCTATCCT/bigwig/all_reads.bw'.format(merge=merge) for merge in params_list],
        #bw_merge_H3K27me3 =['{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads.bw'.format(merge=merge) for merge in params_list],
        #bam_merge_ATAC =['{merge}/ATAC_TATAGCCT/cellranger/outs/possorted_bam.bam'.format(merge=merge) for merge in params_list],
        #bam_merge_H3K27ac =['{merge}/H3K27ac_CCTATCCT/cellranger/outs/possorted_bam.bam'.format(merge=merge) for merge in params_list],
        #bam_merge_H3K27me3 =['{merge}/H3K27me3_ATAGAGGC/cellranger/outs/possorted_bam.bam'.format(merge=merge) for merge in params_list]


ruleorder: run_cellranger_ATAC > run_cellranger_RNA > bam_to_bw > run_macs_broad > barcode_metrics_peaks > barcode_metrics_all > cell_selection 
#merge_bw_ATAC > merge_bw_H3K27ac > merge_bw_H3K27me3 > merge_bam_ATAC > merge_bam_H3K27ac > merge_bam_H3K27me3

#rule demultiplex:
    #input:
        #script=workflow.basedir + '/scripts/debarcode.py',
        #fastq=lambda wildcards: glob.glob(config['samples'][wildcards.sample]['fastq_path'] + '/**/*{lane}*R[123]*.fastq.gz'.format(lane=wildcards.lane),recursive=True)
    #output:
        #'{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R1_{suffix}',
        #'{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R2_{suffix}',
        #'{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/{sample}_{number}_{lane}_R3_{suffix}',
    #params:
        #nbarcodes=lambda wildcards: len(config['samples'][wildcards.sample]['barcodes']),
        #out_folder=lambda wildcards: '{sample}/{modality}_{barcode}/fastq/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode),
    #conda: '../envs/nanoscope_debarcode.yaml'
    #shell:
        #"python3 {input.script} -i {input.fastq} -o {params.out_folder} --single_cell --barcode {wildcards.barcode} -p AGACGCTGCCGACGACAGACGCG --name {wildcards.sample} 2>&1"

rule run_cellranger_ATAC:
    input:
        lambda wildcards: get_fastq_for_cellranger(config['samples'][wildcards.sample]['fastq_path'],sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    output:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        frag='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz',
        meta='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        peaks='{sample}/{modality}_{barcode}/cellranger/outs/peaks.bed',
    params:
        cellranger_software=config['general']['cellranger_software'],
        cellranger_ref=config['general']['cellranger_ref'],
        fastq_folder=lambda wildcards: os.getcwd() + '/{sample}/{modality}_{barcode}/fastq/barcode_{barcode}/'.format(sample=wildcards.sample,modality=wildcards.modality,barcode=wildcards.barcode)
    threads: 20
    resources:
        mem_mb = 32000,
        mem_gb=32
    shell:
        'rm -rf {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/cellranger/; '
        'cd {wildcards.sample}/{wildcards.modality}_{wildcards.barcode}/; '
        '{params.cellranger_software} count --id cellranger --reference {params.cellranger_ref} --chemistry=ARC-v1 --fastqs {params.fastq_folder} --localcores={threads}  --localmem={resources.mem_gb}'

rule run_cellranger_RNA:
    input:
        lambda wildcards: get_fastq_for_cellranger_rna(config['samples'][wildcards.sample]['fastq_path_RNA']+ '/**/*{lane}*R[12]*.fastq.gz',sample=wildcards.sample)
    output:
        bam_RNA='{sample}/RNA_AAAAGGGG_H3K27ac/cellranger/outs/possorted_genome_bam.bam'
    params:
        cellranger_software=config['general']['cellranger_software_RNA'],
        cellranger_ref=config['general']['cellranger_ref_RNA'],
        fastq_folder=lambda wildcards: config['samples'][wildcards.sample]['fastq_path_RNA']#+'/02-FASTQ/20231116_LH00217_0027_A22FMT5LT3/'
    threads: 20
    resources:
        mem_mb = 32000
    shell:
        'rm -rf {wildcards.sample}/RNA_AAAAGGGG_H3K27ac/cellranger/; '
        'cd {wildcards.sample}/RNA_AAAAGGGG_H3K27ac/; '
        '{params.cellranger_software} count --id cellranger --transcriptome {params.cellranger_ref} --chemistry=ARC-v1 --fastqs {params.fastq_folder}'

rule bam_to_bw: # For QC reasons
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        bigwig='{sample}/{modality}_{barcode}/bigwig/all_reads.bw'
    threads: 20
    conda: config['general']['workflow_dir'] + 'envs/nanoscope_deeptools.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule run_macs_broad:
    input:
        cellranger_bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam'
    output:
        broad_peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'
    params:
        macs_outdir='{sample}/{modality}_{barcode}/peaks/macs_broad/',
        macs_genome=config['general']['macs_genome']
    conda: config['general']['workflow_dir'] +'envs/nanoscope_deeptools.yaml'
    resources:
        mem_mb = 16000
    shell:
        'macs2 callpeak -t {input} -g {params.macs_genome} -f BAMPE -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup 1 --broad-cutoff 0.1 '
        '--max-gap 1000 --broad 2>&1 '

rule barcode_metrics_peaks:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
    output:
        overlap='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
    params:
        get_cell_barcode=config['general']['workflow_dir'] + 'workflow//scripts/get_cell_barcode.awk',
        add_sample_to_list=config['general']['workflow_dir'] + 'workflow//scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    conda: config['general']['workflow_dir'] + 'envs/nanoscope_deeptools.yaml'
    resources:
        mem_mb = 16000
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
        'sort -T {params.tmpdir} | uniq -c > {output.overlap} && [[ -s {output.overlap} ]] ; '

rule barcode_metrics_all:
    input:
        bam='{sample}/{modality}_{barcode}/cellranger/outs/possorted_bam.bam',
    output:
        all_bcd='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt'
    params:
        get_cell_barcode=config['general']['workflow_dir'] + 'workflow//scripts/get_cell_barcode.awk',
        add_sample_to_list=config['general']['workflow_dir']+ 'workflow//scripts/add_sample_to_list.py',
        tmpdir=config['general']['tempdir']
    conda: config['general']['workflow_dir'] + 'envs/nanoscope_deeptools.yaml'
    resources:
        mem_mb = 16000
    shell:
        'mkdir -p {params.tmpdir}; '
        ' samtools view -f2 {input.bam}| '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" |  '
        'sort -T {params.tmpdir} | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] ; '

rule cell_selection:
    input:
        bcd_all='{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt',
        bcd_peak='{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt',
        peaks='{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak',
        metadata='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv',
        fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'
    output:
        '{sample}/{modality}_{barcode}/cell_picking/cells_10x.png',
        '{sample}/{modality}_{barcode}/cell_picking/cells_picked.png',
        '{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
    params:
        script=config['general']['workflow_dir'] + 'workflow/scripts/pick_cells.R',
        out_prefix='{sample}/{modality}_{barcode}/cell_picking/',
    resources:
        mem_mb = 25000
    conda: config['general']['workflow_dir'] + 'envs/nanoscope_pick_cells.yaml'
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --modality {wildcards.modality} --min_reads 0.5  --sample {wildcards.sample} --out_prefix {params.out_prefix}"


#rule generate_csv:
    #input:
        #script=workflow.basedir + '/scripts/csv_generator.py',
        #metadata='{sample}/{modality}_{barcode}/cellranger/outs/singlecell.csv'.format(sample=sample,modality=modality,barcode=barcode),
        #fragments='{sample}/{modality}_{barcode}/cellranger/outs/fragments.tsv.gz'.format(sample=sample,modality=modality,barcode=barcode)
    #output:
        #csv="/date/gcb/gcb_MZ/multiNanoCT/merge/{modality}_{barcode}/library.csv"
    #params:
        #sample_data=lambda wildcards:config['samples'][wildcards.sample]
    #conda: '../envs/nanoscope_debarcode.yaml'
    #shell:
        #"python3 {input.script} -o {output.csv} --sample_data {wildcards.sample} --fragments_path {input.fragments} --cells_path {input.metadata} 2>&1"

#rule aggregate_fragments:
    #input:
        #csv="/date/gcb/gcb_MZ/multiNanoCT/samples/FC_Droplet_Paired-Tag/{merge}/{modality_n}_{barcode_n}/library.csv" 
    #output:
        #fragments="{merge}/{modality_n}_{barcode_n}/cellranger/outs/fragments.tsv.gz",
        #peaks="{merge}/{modality_n}_{barcode_n}/cellranger/outs/peaks.bed"
    #params:
       # cellranger_software=config['general']['cellranger_software'],
        #cellranger_ref=config['general']['cellranger_ref'],
       # normalization=config['general']['norm']
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #'rm -rf {wildcards.merge}/{wildcards.modality_n}_{wildcards.barcode_n}/cellranger/; '
        #'cd {wildcards.merge}/{wildcards.modality_n}_{wildcards.barcode_n}/; '
        #'{params.cellranger_software} aggr --id=cellranger --reference={params.cellranger_ref} --normalize={params.normalization} --csv={input.csv} --localcores={threads}  --localmem={resources.mem_gb}'

#rule merge_bw_ATAC:
    #input:
        #bw=['{sample}/ATAC_TATAGCCT/bigwig/all_reads.bw'.format(sample=sample) for sample in samples_list]
    #output:
        #bdg ='{merge}/ATAC_TATAGCCT/bigwig/all_reads.bdg',
        #sorted = '{merge}/ATAC_TATAGCCT/bigwig/all_reads_sorted.bdg',
        #bw_merge = "{merge}/ATAC_TATAGCCT/bigwig/all_reads.bw"
    #params:
        #chrom_sizes = config['general']['chrom_sizes']
    #conda: '/home/mattia/miniconda3/envs/bedtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #/home/mattia/UCSC/bigWigMerge {input.bw} {output.bdg}
        
        #bedtools sort -i {output.bdg} > {output.sorted}
        
        #/home/mattia/UCSC/bedGraphToBigWig {output.sorted} {params.chrom_sizes} {output.bw_merge}
        #"""
#rule merge_bw_H3K27ac:
    #input:
        #bw=['{sample}/H3K27ac_CCTATCCT/bigwig/all_reads.bw'.format(sample=sample) for sample in samples_list]
    #output:
        #bdg ='{merge}/H3K27ac_CCTATCCT/bigwig/all_reads.bdg',
        #sorted = '{merge}/H3K27ac_CCTATCCT/bigwig/all_reads_sorted.bdg',
        #bw_merge = "{merge}/H3K27ac_CCTATCCT/bigwig/all_reads.bw"
    #params:
        #chrom_sizes = config['general']['chrom_sizes']
    #conda: '/home/mattia/miniconda3/envs/bedtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #/home/mattia/UCSC/bigWigMerge {input.bw} {output.bdg}
        
        #bedtools sort -i {output.bdg} > {output.sorted}
        
        #/home/mattia/UCSC/bedGraphToBigWig {output.sorted} {params.chrom_sizes} {output.bw_merge}
        #"""
#rule merge_bw_H3K27me3:
    #input:
        #bw=['{sample}/H3K27me3_ATAGAGGC/bigwig/all_reads.bw'.format(sample=sample) for sample in samples_list]
    #output:
        #bdg ='{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads.bdg',
        #sorted = '{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads_sorted.bdg',
        #bw_merge = "{merge}/H3K27me3_ATAGAGGC/bigwig/all_reads.bw"
    #params:
        #chrom_sizes = config['general']['chrom_sizes']
    #conda: '/home/mattia/miniconda3/envs/bedtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #/home/mattia/UCSC/bigWigMerge {input.bw} {output.bdg}
        
        #bedtools sort -i {output.bdg} > {output.sorted}
        
        #/home/mattia/UCSC/bedGraphToBigWig {output.sorted} {params.chrom_sizes} {output.bw_merge}
        #"""
#rule merge_bam_ATAC:
    #input:
        #bam=['{sample}/ATAC_TATAGCCT/cellranger/outs/possorted_bam.bam'.format(sample=sample) for sample in samples_list]
    #output:
        #merged ='{merge}/ATAC_TATAGCCT/cellranger/outs/possorted_bam.bam'
    #conda: '/home/mattia/miniconda3/envs/samtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #samtools merge {output.merged} {input.bam} -@ {threads}

        #samtools index {output.merged} -@ {threads}

        #"""
#rule merge_bam_H3K27ac:
    #input:
        #bam=['{sample}/H3K27ac_CCTATCCT/cellranger/outs/possorted_bam.bam'.format(sample=sample) for sample in samples_list]
    #output:
        #merged = "{merge}/H3K27ac_CCTATCCT/cellranger/outs/possorted_bam.bam"
    #conda: '/home/mattia/miniconda3/envs/samtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
        #"""
        
        #samtools merge {output.merged} {input.bam}  -@ {threads}

        #samtools index {output.merged} -@ {threads}

        #"""
#rule merge_bam_H3K27me3:
    #input:
        #bam=['{sample}/H3K27me3_ATAGAGGC/cellranger/outs/possorted_bam.bam'.format(sample=sample) for sample in samples_list]
    #output:
        #merged ='{merge}/H3K27me3_ATAGAGGC/cellranger/outs/possorted_bam.bam'
    #conda: '/home/mattia/miniconda3/envs/samtools.yml'
    #threads: 20
    #resources:
        #mem_mb = 32000,
        #mem_gb = 32
    #shell:
       #"""
        
        #samtools merge {output.merged} {input.bam}  -@ {threads}

        #samtools index {output.merged} -@ {threads}
        #"""