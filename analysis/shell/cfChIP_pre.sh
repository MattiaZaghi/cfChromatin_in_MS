#!/bin/bash

# Define an array of sample names
samples=("07068-P-MS-Rituximab-Stable_H3K27ac_ChIP" "H9-P-Ctrl_H3K27ac_ChIP"
"12179-P-MS-Rituximab-Progressive_H3K27ac_ChIP" "P12179-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP"
"14131-P-MS-Rituximab-Stable_H3K27ac_ChIP" "P14020-P-MS-Rituximab-Stable_H3K27ac_ChIP"
"18070-P-MS-Rituximab-Progressive_H3K27ac_ChIP" "P14131-P-MS-Rituximab-Stable-pA_H3K27ac_ChIP"
"18075-P-MS-Rituximab-Progressive_H3K27ac_ChIP" "P14229-P-MS-Rituximab-Stable-pA_H3K27ac_ChIP"
"H10-P-Ctrl_H3K27ac_ChIP" "P14245-P-MS-Rituximab-Stable_H3K27ac_ChIP"
"H11-P-Ctrl_H3K27ac_ChIP" "P15024-P-MS-Rituximab-Stable_H3K27ac_ChIP"
"H12-P-Ctrl_H3K27ac_ChIP" "P15041-P-MS-Rituximab-Stable_H3K27ac_ChIP"
"H13-P-Ctrl_H3K27ac_ChIP" "P16186-P-MS-Rituximab-Stable_H3K27ac_ChIP"
"H14-P-Ctrl_H3K27ac_ChIP" "P18070-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP"
"H15-P-Ctrl_H3K27ac_ChIP" "P18075-P-MS-Rituximab-Prog-pA_H3K27ac_ChIP"
"H16-P-Ctrl_H3K27ac_ChIP" "P20015-P-Ctrl_H3K27ac_ChIP"
"H17-P-Ctrl_H3K27ac_ChIP" "P20027-P-Ctrl_H3K27ac_ChIP"
"H18-P-Ctrl_H3K27ac_ChIP" "P20030-P-Ctrl_H3K27ac_ChIP"
"H19-P-Ctrl_H3K27ac_ChIP" "P20040-P-Ctrl_H3K27ac_ChIP"
"H20-P-Ctrl_H3K27ac_ChIP" "P24116-P-MS-New-RR_H3K27ac_ChIP"
"H21-P-Ctrl_H3K27ac_ChIP" "P24117-P-MS-New-RR_H3K27ac_ChIP"
"H22-P-Ctrl_H3K27ac_ChIP" "P24118-P-MS-New-RR_H3K27ac_ChIP"
"H23-P-Ctrl_H3K27ac_ChIP" "P24126-P-MS-New-RR_H3K27ac_ChIP"
"H24-P-Ctrl_H3K27ac_ChIP" "P24132-P-MS-New-RR_H3K27ac_ChIP"
"H5-P-Ctrl_H3K27ac_ChIP" "P24134-P-MS-New-RR_H3K27ac_ChIP"
"H6-P-Ctrl_H3K27ac_ChIP" "P24136-P-MS-New-RR_H3K27ac_ChIP"
"H7-P-Ctrl_H3K27ac_ChIP")

# Export the necessary variables for xargs
export Analysis=Analysis
export Output=Analysis/Output/H3K36me3
export m=H3K36me3
export track=Analysis/Tracks/H3K36me3

# Define the function to run the R script
run_script() {
    sample=$1
    # First step create rdata computing counts and background
    Rscript --vanilla ${Analysis}/cfChIP-seq/ProcessBEDFiles_MZ_2.R  ${sample}
    # Additional step write meta plot, normcounts table and backgroundplots
}

# Export the function
export -f run_script

# Use xargs to run the function in parallel
printf '%s\n' "${samples[@]}" | xargs -n 50 -P 50 -I {} bash -c 'run_script "$@"' _ {}


