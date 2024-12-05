#!/usr/bin/env python3

import json
import os
import csv
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", help="Required. the FULL path to the fastq folder")
parser.add_argument("--meta", help="Required. the FULL path to the tab delimited meta file")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"
assert args.meta is not None, "please provide the path to the meta file"

# Collect all the fastq.gz full paths into a list
fastq_paths = []

for root, dirs, files in os.walk(args.fastq_dir):
    for file in files:
        if file.endswith(".gz"):
            full_path = join(root, file)
            fastq_paths.append(full_path)

FILES = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

with open(args.meta, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    # Skip the header
    header = next(reader)
    for row in reader:
        print(f"Processing row: {row}")  # Debugging line
        if len(row) < 6:
            print(f"Skipping incomplete row: {row}")
            continue
        sample_name = row[0].strip()
        fastq_name = row[1].strip()
        factor = row[2].strip()
        assay = row[3].strip()
        sample_type = row[4].strip()
        reference = row[5].strip()
        # Assume the file name in the metafile is contained in the fastq file path
        fastq_full_path = [x for x in fastq_paths if fastq_name in x]
        if fastq_full_path:
            FILES[sample_name][factor][assay].extend(fastq_full_path)
        else:
            print(f"Sample {sample_name} missing {sample_type} {fastq_name} fastq files")

print()
sample_num = len(FILES.keys())
print(f"Total {sample_num} unique samples will be processed")
print("------------------------------------------")
for sample_name in sorted(FILES.keys()):
    for factor in FILES[sample_name]:
        for assay in FILES[sample_name][factor]:
            fastq_file = "".join(FILES[sample_name][factor][assay])
            print(f"Sample {sample_name}'s {factor} fastq path is {fastq_file}")
print("------------------------------------------")
for sample in FILES.keys():
    print(f"{sample} has {len(FILES[sample])} marks")
print("------------------------------------------")
print("Check the samples.json file for fastqs belonging to each sample")
print()
js = json.dumps(FILES, indent=4, sort_keys=True)
with open('samples_dag.json', 'w') as outfile:
    outfile.write(js)
