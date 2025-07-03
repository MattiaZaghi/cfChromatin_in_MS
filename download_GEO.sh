#!/usr/bin/env bash
###############################################################################
# get_cfchip_GSE243474.sh
# Download, extract and rename cfChIP-seq H3K4me3 / H3K27ac BED files from
# GEO series GSE243474, while printing progress messages.
#
# Usage:
#   bash get_cfchip_GSE243474.sh        # normal, concise logging
#   bash get_cfchip_GSE243474.sh -v     # verbose (all commands echoed)
###############################################################################
set -euo pipefail

###############################################################################
# 0.  Logging helpers
###############################################################################
log() { printf '[%(%F %T)T] %s\n' -1 "$*"; }          # date-stamped message
[[ ${1:-} == "-v" ]] && { log "Verbose mode on"; set -x; }

###############################################################################
# 1.  Constants
###############################################################################
ACC=GSE243474
SERIES_DIR=${ACC:0:6}nnn                 #  GSE243nnn  (GEO FTP convention)
FTP=ftp://ftp.ncbi.nlm.nih.gov/geo/series/${SERIES_DIR}/${ACC}
RAW=${ACC}_RAW.tar
OUTDIR=${ACC}_cfChIP

###############################################################################
# 2.  Prepare working directory
###############################################################################
mkdir -p "$OUTDIR"
cd "$OUTDIR"

###############################################################################
# 3.  Download RAW tarball
###############################################################################
log "Step 1/4 – Downloading ${RAW} from GEO …"
if [[ -f $RAW ]]; then
    log "  • ${RAW} already present – skipping download."
else
    wget -nv -c "${FTP}/suppl/${RAW}"
fi

###############################################################################
# 4.  Select BED files of interest (only *.bed.gz, no *_sorted_peaks*)
###############################################################################
log "Step 2/4 – Scanning archive for H3K4me3 / H3K27ac cfChIP BEDs …"
tar -tf "$RAW" \
 | grep -E '(_K4[^/]*\.bed\.gz|_K27[^/]*\.bed\.gz)$' \
 | grep -v 'narrowPeak' \
 > wanted.txt
   # keeps only pure BEDs
NUM_WANTED=$(wc -l < wanted.txt)
[[ $NUM_WANTED -eq 0 ]] && { log "  • No matching files – aborting."; exit 1; }
log "  • Found ${NUM_WANTED} matching files."

###############################################################################
# 5.  Extract only those files  (flatten every path component)
###############################################################################
log "Step 3/4 – Extracting selected files …"
tar -xvf "$RAW" -T wanted.txt --transform='s:.*/::'
log "  • Extraction completed."



###############################################################################
# 6.  Rename files: CTRL = healthy, TUMOR = patient + tumour-type + counter
###############################################################################
log "Step 4/4 – Renaming files …"                                          # [1]
shopt -s nullglob                                                          # [2]

declare -A seen                                                            # [2]
for f in *_K4*.bed.gz *_K27*.bed.gz; do                                    # [1]
    id=${f%%_*}                    # GSM7789068, HP01, mLC33 …            # [1]
    id=${id#TUMOR_}                # drop any stray leading word           # [2]
    mark=$(grep -oE 'K4[^_]*|K27[^_]*' <<<"$f")                            # [1]

    # NEW – grab the tumour-type (token right after the histone mark)      # [1]
    tumour_type=$(echo "$f" | awk -F'_' -v m="$mark" '{                    \
                    for(i=1;i<=NF;i++) if($i==m){print $(i+1); break}      \
                  }')                                                      # [2]

    if [[ $id =~ ^HP ]]; then                                             # control sample  # [1]
        new="CTRL_${id}_${mark}.bed.gz"                                    # [1]
    else                                                                   # tumour sample   # [1]
        key="${tumour_type}_${mark}"                                       # [2]
        seen[$key]=$(( ${seen[$key]:-0} + 1 ))                             # [2]
        idx=$(printf "R%02d" "${seen[$key]}")                              # [2]
        mark=${mark/_ac/}                 # shorten K27ac → K27            # [1]
        new="TUMOR_${tumour_type}_${mark}_${idx}.bed.gz"                   # [1]
    fi
    mv -v -- "$f" "$new"                                                   # [2]
done

log "  • Finished renaming – $(ls CTRL_*.bed.gz | wc -l) controls and \
$(ls TUMOR_*.bed.gz | wc -l) tumours."                                     # [1]



###############################################################################
# 7.  Wrap-up
###############################################################################
log "All done – result files are in: $(pwd)"
