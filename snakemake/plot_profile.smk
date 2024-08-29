# ------ sample handling ----------
BED="TSS_cfChromatin".split() # Regions on which to perform heatmaps and profiles plot
CHIP="Plasma_H3K27ac".split() # Bigwig tracks to plot H3K27ac_CCTATCCT 
STAT="median".split()# Value to plot
ZMAX= "Auto".split() # Max Value to plot in heatmaps
KMEANS="".split() # Number of Kmeans to plot into heatmaps 


rule all:
    input:
        expand('/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{chip}_{bed}_{stat}_Compute_Matrix_profile.png',bed=BED,chip=CHIP,stat=STAT),
        expand('/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{chip}_{bed}_{zmax}_{stat}_Compute_Matrix_heatmap.png', chip=CHIP,bed=BED,stat=STAT, zmax=ZMAX)
        
rule Compute_Matrix:
    input:
        bigwig1='/date/gcb/gcb_MZ/Round2_fresh_hg19/coverage/deeptools/R/H1-P_H3K27ac_ChIP.bw',
        bigwig2='/date/gcb/gcb_MZ/Round2_fresh_hg19/coverage/deeptools/R/H2-P_H3K27ac_ChIP.bw',
        bigwig3='/date/gcb/gcb_MZ/Round2_fresh_hg19/coverage/deeptools/R/H3-P_H3K27ac_ChIP.bw',
        bigwig4='/date/gcb/gcb_MZ/Round2_fresh_hg19/coverage/deeptools/R/H4-P_H3K27ac_ChIP.bw',
        BED= '/date/gcb/gcb_MZ/Analysis/cfChIP-seq/SetupFiles/H3K27ac/TSS_cfChromatin.bed'
    conda: "/home/mattia/miniconda3/envs/deeptools.yml"
    resources:
        mem_mb = 16000
    threads: 5
    output:
        Compute_Matrix='/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{chip}_{bed}_{stat}_Compute_Matrix'
    shell:
        "computeMatrix reference-point -S {input.bigwig1}  {input.bigwig2} {input.bigwig3} -R {input.BED} -o {output.Compute_Matrix} -a 2000 -b 7000 -p {threads} --averageTypeBins {wildcards.stat} --referencePoint center"
        

rule Plot_profile_automax:
    input:
        Compute_Matrix='/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{bed}_{stat}_Compute_Matrix'
    output:
        profile_automax="/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{bed}_{stat}_Compute_Matrix_profile.png"
    conda: "/home/mattia/miniconda3/envs/deeptools.yml"
    resources:
        mem_mb = 16000
    threads: 5
    shell:
        """
        plotProfile --matrixFile {input.Compute_Matrix} --outFileName {output.profile_automax} --dpi 300 --averageType {wildcards.stat} --plotWidth 7.5 --perGroup --colors '#FF9900' '#FF6633' '#3399FF' --samplesLabel '' '' '' --startLabel "5'" --endLabel "3'" --regionsLabel '{wildcards.bed}'
        """
        
rule Plot_heatmap_kmeans: #this is the code used to obtain clusters of differential chromatin accessibility in different samples and cell types
    input:
        Compute_Matrix='/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{bed}_{stat}_Compute_Matrix'
    output:
        Heatmap_kmeans="/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{bed}_{zmax}_{stat}_Compute_Matrix_heatmap.png",
        Heatmap_bed= "/date/gcb/gcb_MZ/Round2_fresh_hg19/deeptools/{bed}_{zmax}_{stat}_Compute_Matrix_heatmap.bed"
    conda: "/home/mattia/miniconda3/envs/deeptools.yml"
    resources:
        mem_mb = 16000
    threads: 5
    shell:
        """
        plotHeatmap --matrixFile {input.Compute_Matrix} --outFileName {output.Heatmap_kmeans} --heatmapWidth 3.5 --heatmapHeight 20 --dpi 300 --sortRegions descend --averageTypeSummaryPlot {wildcards.stat} --refPointLabel center --colorMap Greys  --missingDataColor  lightgrey --whatToShow  'heatmap and colorbar' --samplesLabel "H1-P" "H2-P" "H3-P"  --startLabel "5'" --endLabel "3'" --legendLocation none --outFileSortedRegions  {output.Heatmap_bed}  --zMax {wildcards.zmax} 
        """