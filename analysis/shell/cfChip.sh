#!/bin/bash

# Define an array of sample names
samples=("H1-P_H3K36me3_ChIP" "H2-P_H3K36me3_ChIP" "H3-P_H3K36me3_ChIP" "H4-P_H3K36me3_ChIP" "H003.6" "H012.1" "H012.2" "H013.1" "H013.2")

# Export the necessary variables for xargs
export Analysis=Analysis
export Output=Analysis/Output/H3K36me3
export m=H3K36me3
export track=Analysis/Tracks/H3K36me3

# Define the function to run the R script
run_script() {
    sample=$1
    # First step create rdata computing counts and background
    Rscript --vanilla ${Analysis}/cfChIP-seq/ProcessBEDFiles.R -r ${Analysis} --outputdir ${Output} --trackdir ${track} -m ${m} -BCN ${sample}
    # Additional step write meta plot, normcounts table and backgroundplots
    Rscript --vanilla ${Analysis}/cfChIP-seq/ProcessBEDFiles.R -r ${Analysis} --outputdir ${Output} --trackdir ${track} -m ${m} --backgroundplot --normcounts=NORMCOUNTS_${sample} ${sample}
    # Additional step write tracks
    Rscript --vanilla ${Analysis}/cfChIP-seq/ProcessBEDFiles.R -r ${Analysis} --outputdir ${Output} --trackdir ${track} -m ${m} -T ${sample}
}

# Export the function
export -f run_script

# Use xargs to run the function in parallel
printf '%s\n' "${samples[@]}" | xargs -n 40 -P 40 -I {} bash -c 'run_script "$@"' _ {}


# Export the function
export -f run_script

# Use xargs to run the function in parallel
printf '%s\n' "${samples[@]}" | xargs -n 40 -P 40 -I {} bash -c 'run_script "$@"' _ {}
