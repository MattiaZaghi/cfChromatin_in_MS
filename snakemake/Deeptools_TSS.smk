
# ------ sample handling ----------
BED="House_keeping".split() # Regions on which to perform heatmaps and profiles plot
CHIP="cfH3K27ac_all".split() # Bigwig tracks to plot
STAT="median".split()# Value to plot
ZMAX= "auto".split() # Max Value to plot in heatmaps
KMEANS="1".split() # Number of Kmeans to plot into heatmaps 


rule all:
    input:
        #expand('{chip}_in_{bed}_{stat}_Compute_Matrix_profile.png',bed=BED,chip=CHIP,stat=STAT),
        expand('{chip}_in_{bed}_{zmax}_{stat}_Compute_Matrix_{kmeans}_heatmap.png', bed=BED,chip=CHIP,stat=STAT, kmeans=KMEANS,zmax=ZMAX)
        
rule Compute_Matrix:
    input:
        bigwig1='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/GSM7787985_HP038748_H3K27Ac.bw',
        bigwig2='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/GSM7788108_K27_HPC81_2.bw',
        bigwig3='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/H23-P-Ctrl_H3K27ac_ChIP.bw',
        bigwig4='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/P18070-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP.bw',
        bigwig5='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP.bw',
        bigwig6='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/16057-P-100_H3K27ac_ChIP.bw',
        bigwig7='/proj/user/mattia/Analysis/Tracks/H3K27ac_hg38/16170-C_H3K27ac_ChIP.bw',
        BED= '/proj/user/mattia/Analysis/cfChIP-seq/SetupFiles/H3K27ac_hg38/{bed}.bed'
    output:
        Compute_Matrix='{chip}_in_{bed}_{stat}_Compute_Matrix'
    resources:
        mem_mb = 10000
    threads: 30
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        "computeMatrix  reference-point -S {input.bigwig1} {input.bigwig2}  {input.bigwig3} {input.bigwig4} {input.bigwig5} {input.bigwig6} {input.bigwig7} -R {input.BED} -o {output.Compute_Matrix} -a 3000  -b 3000 -p {threads} --averageTypeBins {wildcards.stat} "
        

rule Plot_profile_automax:
    input:
        Compute_Matrix='{chip}_in_{bed}_{stat}_Compute_Matrix'
    output:
        profile_automax='{chip}_in_{bed}_{stat}_Compute_Matrix_profile.png'
    resources:
        mem_mb = 10000
    threads: 1
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        plotProfile --matrixFile {input.Compute_Matrix} --outFileName {output.profile_automax} --dpi 300 --averageType {wildcards.stat} --plotWidth 20 --perGroup --colors  '#fcba03' '#3d85c6'  --samplesLabel '' ''  --startLabel "TSS" --endLabel "TES" --regionsLabel '{wildcards.chip} in {wildcards.bed}'
        """
        # ChIP=#d33928' S3='#fcba03' NanoCT'#3d85c6'
        #'#1f77b4' '#ff7f0e' '#2ca02c' '#9467bd' '#8c564b'
rule Plot_heatmap_kmeans: #this is the code used to obtain clusters of differential chromatin accessibility in different samples and cell types
    input:
        Compute_Matrix='{chip}_in_{bed}_{stat}_Compute_Matrix'
    output:
        Heatmap_kmeans="{chip}_in_{bed}_{zmax}_{stat}_Compute_Matrix_{kmeans}_heatmap.png",
        Heatmap_bed= "{chip}_in_{bed}_{zmax}_{stat}_Compute_Matrix_{kmeans}_heatmap.bed"
    resources:
        mem_mb = 10000
    threads: 1
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        plotHeatmap --matrixFile {input.Compute_Matrix} --outFileName {output.Heatmap_kmeans} --heatmapWidth 8 --heatmapHeight 20 --dpi 300 --sortRegions no --averageTypeSummaryPlot {wildcards.stat} --kmeans {wildcards.kmeans}  --refPointLabel center --colorMap Oranges Oranges Oranges Oranges Oranges Oranges Oranges  --missingDataColor  lightgrey --whatToShow  'plot, heatmap and colorbar' --samplesLabel  "Baca-Ctrl" "Baca-Canc" "Plasma V3 Fresh" "Plasma V3" "Plasma V2" "Plasma V1" "CSF V2"  --startLabel "TSS" --endLabel "TES" --legendLocation none --outFileSortedRegions  {output.Heatmap_bed}  --zMin 0 0 0 0 0 0 0 --zMax 4 4 1 1 1 1 1
        """

