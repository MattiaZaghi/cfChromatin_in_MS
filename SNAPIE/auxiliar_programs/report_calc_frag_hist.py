#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2023 @chris-cheshire
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Author: @chris-cheshire
# Modified by: @prc992

import os
import glob
import argparse

import numpy as np
import pandas as pd

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = "Calclate fragment histogram"

parser = argparse.ArgumentParser(description=Description)

## REQUIRED PARAMETERS
parser.add_argument("--frag_path")
parser.add_argument("--output")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

print("Calclate fragment histogram")

# Init
frag_path = os.path.abspath(args.frag_path)
frag_hist = None

# Create list of deeptools raw fragment files
dt_frag_list = glob.glob(frag_path)
dt_frag_list.sort()

for i in list(range(len(dt_frag_list))):
    # Create dataframe from csv file for each file and save to a list of data frames
    dt_frag_i = pd.read_csv(dt_frag_list[i], sep="\t", names=["Size", "Occurrences","Sample"],skiprows=2)
    dt_frag_i = dt_frag_i.drop('Sample', axis=1)
    frag_base_i = os.path.basename(dt_frag_list[i])

    # Split txt file paths on dots
    sample_id_list = frag_base_i.split(".")

    # Join list on the elements of the sample id
    separator = ""
    sample_id = separator.join(sample_id_list[0:-2])

    # Round column and convert occurrences to int
    dt_frag_i = dt_frag_i.round(1)
    dt_frag_i = dt_frag_i.astype(int)

    # Create long forms of fragment histograms
    dt_frag_i_long = np.repeat(dt_frag_i["Size"].values, dt_frag_i["Occurrences"].values)
    dt_sample_i_long = np.repeat(sample_id, len(dt_frag_i_long))
    dt_sample_i_short = np.repeat(sample_id, dt_frag_i.shape[0])

    if i == 0:
        frags_arr = dt_frag_i_long
        sample_arr = dt_sample_i_long
        sample_short = dt_sample_i_short
        frag_hist = dt_frag_i
    else:
        frags_arr = np.append(frags_arr, dt_frag_i_long)
        sample_arr = np.append(sample_arr, dt_sample_i_long)
        sample_short = np.append(sample_short, dt_sample_i_short)
        frag_hist = frag_hist.append(dt_frag_i)

frag_hist["sample"] = sample_short
frag_hist = frag_hist.reset_index(drop=True)

size_list = frag_hist["Size"].to_numpy().astype(str)
occurrences_list = frag_hist["Occurrences"].to_numpy().astype(str)
size_list_sep = np.core.defchararray.add(size_list, " : ")
x_y_list = np.core.defchararray.add(size_list_sep, occurrences_list)

sample_grp = frag_hist[["sample"]].groupby(["sample"]).size().reset_index()
first_line = "data:"

for i in list(range(sample_grp.shape[0])):
    sample_i = sample_grp.at[i, "sample"]
    str_list = x_y_list[(frag_hist["sample"] == sample_i)]

    x_y_str = ", ".join(str_list)
    full_line_i = "    '" + sample_i + "' : {" + x_y_str + "}"
    if i == 0:
        frag_len_hist_mqc_dict = "\n".join([first_line, full_line_i])

    else:
        frag_len_hist_mqc_dict = "\n".join([frag_len_hist_mqc_dict, full_line_i])

txt_mqc = open(args.output, "w")
txt_mqc.write(frag_len_hist_mqc_dict)
txt_mqc.close()