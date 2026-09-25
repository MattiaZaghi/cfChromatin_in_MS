#!/bin/bash
#SBATCH --job-name=pb_markers
#SBATCH --output=pb_markers_%j.out
#SBATCH --error=pb_markers_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
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

BASE=/cfs/klemming/projects/supr/uppstore2017150/anita

python3 /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/build_tissue_markers_sadeh.py \
    --no-gc-correction \
    --window-types TSS,EXTRA_TSS,Enhancer \
    --top-n 400 --high 22.6 --low 15.3 \
    --blood-veto-frac 0.25 --min-markers-floor 20 \
    --top-n-overrides brain=2290,breast=387,digestive=567,heart=1000,prostate=360,t_cell=50,neutrophil=700,lung=250,liver=2000,gi_smooth_muscle=2000,blood_vessel=50,skin=50,skeletal_muscle_f=400,megakaryocyte_mep=100,adipose=400,nk=100,b_cell=100 \
    --out /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/tissue_markers_sadeh

# for brain
# --top-n 400 --high 22.6 --low 15.3
# --top-n-overrides brain=2290,breast=387,digestive=445,heart=20,prostate=360,t_cell=436,neutrophil=700,lung=50,liver=50,t_cell=50,gi_smooth_muscle=100,blood_vessel=50,skin=50,skeletal_muscle_f=50,gi_smooth_muscle=50,digestive=100,megakaryocyte_mep=100,adipose=100,nk=100,b_cell=100
# for brain_pb
# --top-n 300 --high 22.6 --low 15.3
# --top-n-overrides brain_pb=2000,breast=300,digestive=700,heart=20,prostate=500,t_cell=100,neutrophil=700,lung=50,liver=50,t_cell=50,gi_smooth_muscle=100,blood_vessel=50,skin=50,skeletal_muscle_f=50,gi_smooth_muscle=50,megakaryocyte_mep=100,adipose=100,nk=100,b_cell=100

python3 ${BASE}/cfChromatin_in_MS/tissue_markers/universe_and_anchors.py \
    --markers_dir       ${BASE}/cfChromatin_in_MS/tissue_markers_sadeh \
    --peaks_dir         ${BASE}/cfChromatin_in_MS/MS-cfEpiMap/reference/rre/raw \
    --housekeeping_seed ${BASE}/cfChromatin_in_MS/MS-cfEpiMap/reference/rre/housekeeping_seed.bed \
    --min_anchor_types  0.8 \
    --out ${BASE}/cfChromatin_in_MS/tissue_markers_sadeh

#cp /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/tissue_markers_sadeh_brain_pb/universe.bed /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/tissue_markers_sadeh_brain_pb/universe_top.bed
cp /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/tissue_markers_sadeh/universe.bed /cfs/klemming/projects/supr/uppstore2017150/anita/cfChromatin_in_MS/tissue_markers_sadeh/universe_top.bed