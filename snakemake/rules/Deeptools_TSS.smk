
# ------ sample handling ----------
BED="meta_genes_high".split() # Regions on which to perform heatmaps and profiles plot
CHIP="TSS_comp".split() # Bigwig tracks to plot
STAT="median".split()# Value to plot
ZMAX= "auto".split() # Max Value to plot in heatmaps
KMEANS="1".split() # Number of Kmeans to plot into heatmaps 


rule all:
    input:
        expand('{chip}_in_{bed}_{stat}_Compute_Matrix_profile.png',bed=BED,chip=CHIP,stat=STAT),
        expand('{chip}_in_{bed}_{zmax}_{stat}_Compute_Matrix_{kmeans}_heatmap.png', bed=BED,chip=CHIP,stat=STAT, kmeans=KMEANS,zmax=ZMAX)
        
rule Compute_Matrix:
    input:
        bigwig1='/date/gcb/gcb_MZ/cfChrom_Tak/coverage/deeptools/H32-Ctrl_H3K27ac_ChIP-Tak_RPKM.bw',
        bigwig2='/date/gcb/gcb_MZ/cfChrom_Tak/coverage/deeptools/H44-Ctrl_H3K27ac_ChIP-Tak_RPKM.bw',
        bigwig3='/date/gcb/gcb_MZ/cfChrom_Tak/coverage/deeptools/H43-Ctrl_H3K27ac_ChIP-Tak_RPKM.bw',
        bigwig4='/date/gcb/gcb_MZ/cfChrom_Tak/coverage/deeptools/H30-Ctrl_H3K27ac_ChIP-Tak_RPKM.bw',
        bigwig5='/date/gcb/gcb_MZ/cfChrom_hg19/coverage/deeptools/H24-P-Ctrl_H3K27ac_ChIP-V3_RPKM.bw',
        bigwig6='/date/gcb/gcb_MZ/cfChrom_hg19/coverage/deeptools/H17-P-Ctrl_H3K27ac_ChIP-V3_RPKM.bw',
        bigwig7='/date/gcb/gcb_MZ/cfChrom_hg19/coverage/deeptools/12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP-V2_RPKM.bw',
        bigwig8='/date/gcb/gcb_MZ/cfChrom_hg19/coverage/deeptools/18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP-V2_RPKM.bw',
        BED= '/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/{bed}.bed'
    output:
        Compute_Matrix='{chip}_in_{bed}_{stat}_Compute_Matrix'
    resources:
        mem_mb = 10000
    threads: 30
    conda:
        "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        "computeMatrix  reference-point -S {input.bigwig1} {input.bigwig2} {input.bigwig3} {input.bigwig4} {input.bigwig5} {input.bigwig6} {input.bigwig7} {input.bigwig8} -R {input.BED} -o {output.Compute_Matrix} -a 3000  -b 3000 -p {threads} --averageTypeBins {wildcards.stat} "
        

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
        plotProfile --matrixFile {input.Compute_Matrix} --outFileName {output.profile_automax} --dpi 300 --averageType {wildcards.stat} --plotWidth 20 --perGroup --colors  '#1f77b4' '#ff7f0e' '#2ca02c' '#9467bd'  --samplesLabel '' '' '' '' '' '' --startLabel "TSS" --endLabel "TES" --regionsLabel '{wildcards.chip} in {wildcards.bed}'
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
        plotHeatmap --matrixFile {input.Compute_Matrix} --outFileName {output.Heatmap_kmeans} --heatmapWidth 8 --heatmapHeight 20 --dpi 300 --sortRegions no --averageTypeSummaryPlot {wildcards.stat} --kmeans {wildcards.kmeans}  --refPointLabel center --colorMap Oranges   --missingDataColor  lightgrey --whatToShow  'plot, heatmap and colorbar' --samplesLabel "H32-Tak" "H44-Tak" "H43-Tak" "H30-Tak" "H24-V3" "H17-V3" "12079-V2" "18070-V2" --startLabel "TSS" --endLabel "TES" --legendLocation none --outFileSortedRegions  {output.Heatmap_bed} 
        """

