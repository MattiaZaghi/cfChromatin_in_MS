#config
configfile= "/home/mattia/cfChromatin_in_MS/snakemake/Cut&Tag_bulk/config.yaml"


FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample].keys():
            for factor in FILES[sample].keys():
                MARK_SAMPLES.append(sample + "_" + sample_type+"_"+assay+"_"+factor)

print(MARK_SAMPLES)
