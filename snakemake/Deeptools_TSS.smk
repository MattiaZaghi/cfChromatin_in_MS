
# ------ sample handling ----------
BED="meta_genes_high".split() # Regions on which to perform heatmaps and profiles plot
CHIP="Brain_TSS_samples".split() # Bigwig tracks to plot
STAT="median".split()# Value to plot
ZMAX= "auto".split() # Max Value to plot in heatmaps
KMEANS="1".split() # Number of Kmeans to plot into heatmaps 


rule all:
    input:
        expand('{chip}_in_{bed}_{stat}_Compute_Matrix_profile.png',bed=BED,chip=CHIP,stat=STAT),
        expand('{chip}_in_{bed}_{zmax}_{stat}_Compute_Matrix_{kmeans}_heatmap.png', bed=BED,chip=CHIP,stat=STAT, kmeans=KMEANS,zmax=ZMAX)
        
rule Compute_Matrix:
    input:
        bigwig1='/proj/user/mattia/CPM/H3-P-Ctrl_H3K27ac_ChIP-V2-1D_CPM.bw',
        bigwig2='/proj/user/mattia/CPM/14131-P-MS-Rituximab-Stable_H3K27ac_ChIP-V2_CPM.bw',
        bigwig3='/proj/user/mattia/CPM/12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP-V2_CPM.bw',
        bigwig4='/proj/user/mattia/CPM/18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP-V2_CPM.bw',
        BED= '/proj/user/mattia/CPM/{bed}.bed'
    output:
        Compute_Matrix='{chip}_in_{bed}_{stat}_Compute_Matrix'
    resources:
        mem_mb = 10000
    threads: 30
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        "computeMatrix  reference-point -S {input.bigwig1} {input.bigwig2} {input.bigwig3} {input.bigwig4} -R {input.BED} -o {output.Compute_Matrix} -a 3000  -b 3000 -p {threads} --averageTypeBins {wildcards.stat} "
        

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
        plotProfile --matrixFile {input.Compute_Matrix} --outFileName {output.profile_automax} --dpi 300 --averageType {wildcards.stat} --plotWidth 20 --perGroup --colors  '#1f77b4' '#ff7f0e' '#2ca02c' '#9467bd'  --samplesLabel '' '' '' '' --startLabel "TSS" --endLabel "TES" --regionsLabel '{wildcards.chip} in {wildcards.bed}'
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
        plotHeatmap --matrixFile {input.Compute_Matrix} --outFileName {output.Heatmap_kmeans} --heatmapWidth 8 --heatmapHeight 20 --dpi 300 --sortRegions no --averageTypeSummaryPlot {wildcards.stat} --kmeans {wildcards.kmeans}  --refPointLabel center --colorMap Oranges   --missingDataColor  lightgrey --whatToShow  'plot, heatmap and colorbar' --samplesLabel "Healthy" "MS Stable" "MS progressive 1" "MS progressive 2"  --startLabel "TSS" --endLabel "TES" --legendLocation none --outFileSortedRegions  {output.Heatmap_bed} 
        """

