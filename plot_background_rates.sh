#!/bin/bash
#SBATCH --job-name=plot_bg
#SBATCH --output=plot_bg_%j.out
#SBATCH --error=plot_bg_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --account=naiss2025-22-1791
#SBATCH --partition=shared
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=anita.girardi.1@studenti.unipd.it

# Load conda
module load PDCOLD/21.11
module load Anaconda3/2021.05

source $(conda info --base)/etc/profile.d/conda.sh
# Activate your environment
conda activate /cfs/klemming/projects/supr/uppstore2017150/anita/miniconda3/envs/nanoctar


python3 /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/plot_background_rates.py --fig-width 16 --panel-height 2