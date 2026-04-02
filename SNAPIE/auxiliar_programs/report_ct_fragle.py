#!/usr/bin/env python3

# Initialization
header_file = "ct_header.txt"
data_file = "Fragle.txt"
output_file = "fragle_mqc.csv"

with open(output_file, "w") as out:
    with open(header_file, "r") as hf:
        out.write(hf.read().strip() + "\n")

    with open(data_file, "r") as df:
        for line in df:
            if line.strip():
                out.write(line.strip() + "\n")