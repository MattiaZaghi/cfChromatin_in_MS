#!/usr/bin/env python
# coding: utf-8
import os
import pandas as pd
import argparse

# Path to the header file
header_file = 'signal_header.txt'

# Read the content of the header file
with open(header_file, 'r') as file:
    header_content_original = file.read()


def find_files_with_suffix(root_directory, suffix):
    """
    Find all files in the given directory and its subdirectories that end with the specified suffix.

    :param root_directory: The root directory to search in.
    :param suffix: The suffix to match files against.
    :return: A list of file paths matching the given suffix.
    """
    matching_files = []
    for dirpath, _, filenames in os.walk(root_directory):
        for file in filenames:
            if file.endswith(suffix):
                matching_files.append(os.path.join(dirpath, file))
    return matching_files

current_directory = os.getcwd()

list_files =  find_files_with_suffix(current_directory,'.csv')

data = {
    "File": list_files,
    "Mark": [path.split('/')[-1].split('_Signal_Intensity_Matrix')[0].split('_')[-1] for path in list_files]
}

df = pd.DataFrame(data)

merged_data = {}

for mark, group in df.groupby('Mark'):
    print(f"Processing files for mark: {mark}")

    # Load only non-empty files
    dfs = []
    for file in group['File']:
        if os.path.getsize(file) > 0:  # Ensure file is not empty
            try:
                df_temp = pd.read_csv(file)
                dfs.append(df_temp)
            except pd.errors.EmptyDataError:
                print(f"Warning: {file} is empty. Skipping.")
        else:
            print(f"Warning: {file} is empty. Skipping.")

    if not dfs:  # Skip if no valid data
        print(f"Skipping mark {mark} due to no valid data.")
        continue

    # Concatenate all valid DataFrames
    merged_data[mark] = pd.concat(dfs, ignore_index=True)
    merged_data[mark] = merged_data[mark].drop(columns=['mark'], errors='ignore')

    output_file = f"merged_signal_{mark}_mqc.csv"
    header_content = header_content_original.replace('<MARK>', mark.upper())

    # Write header
    with open(output_file, 'w') as file:
        file.write(header_content + '\n')

    # Write merged data
    merged_data[mark].to_csv(output_file, mode='a', index=False, header=True)