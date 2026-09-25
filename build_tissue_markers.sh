#!/bin/bash
#SBATCH --job-name=encode_bedgraphs
#SBATCH --output=encode_bedgraphs_%j.out
#SBATCH --error=encode_bedgraphs_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --account=naiss2025-22-1791
#SBATCH --partition=shared

# Load conda
module load PDCOLD/21.11
module load Anaconda3/2021.05

source $(conda info --base)/etc/profile.d/conda.sh
# Activate your environment
conda activate /cfs/klemming/projects/supr/uppstore2017150/anita/miniconda3/envs/nanoctar

python3 /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/build_tissue_markers.py --tss-only --top-n 85

