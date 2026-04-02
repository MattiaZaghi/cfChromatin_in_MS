#!/usr/bin/env python
# coding: utf-8
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

# Author: @prc992

import os
import pandas as pd

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

def read_file(str_path):
    try:
        with open(str_path, 'r') as file:
            content = file.read().strip()
    except FileNotFoundError:
        print(f"Error: The file '{header_file}' was not found.")
    except IOError as e:
        print(f"Error: Could not read the file '{header_file}': {e}")
    return content

# Define the suffix
suffix_frags = '_unique_frags.csv'
suffix_peak = '_peaks.narrowPeak'

name_head_file_frags = 'frags_header.txt'
name_head_file_peaks = 'peaks_header.txt'

output_file_frags = "frags_mqc.txt"
output_file_peaks = "peaks_mqc.txt"

current_directory = os.getcwd()

head_file_frags = find_files_with_suffix(current_directory, name_head_file_frags)[0]
head_file_content_frags = read_file(head_file_frags)

head_file_peaks = find_files_with_suffix(current_directory, name_head_file_peaks)[0]
head_file_content_peaks = read_file(head_file_peaks)

frags_files = find_files_with_suffix(current_directory, suffix_frags)    
peak_files = find_files_with_suffix(current_directory, suffix_peak) 

def create_dataframe_from_content_files(file_list,sufix,nameColumn,read_content=True):
    """
    Reads a list of files and returns a DataFrame containing the file name and its content.

    :param file_list: List of file paths to read.
    :return: A pandas DataFrame with columns 'File Name' and 'Content'.
    """
    data = []

    for file_path in file_list:
        try:
            with open(file_path, 'r') as file:
                if read_content == True:
                    content = file.read().strip()
                else:
                    content = sum(1 for _ in file)
                data.append({
                    'SampleName': os.path.basename(file_path).replace(sufix,''),  # Extract only the file name
                    nameColumn: content
                })
        except FileNotFoundError:
            print(f"Error: The file '{file_path}' was not found.")
        except IOError as e:
            print(f"Error: Could not read the file '{file_path}': {e}")

    # Create the DataFrame
    df = pd.DataFrame(data)
    return df

df_frags = create_dataframe_from_content_files(frags_files,suffix_frags,'Fragments',True)
df_peaks = create_dataframe_from_content_files(peak_files,suffix_peak,'Peaks',False)
#merged_df = pd.merge(df_frags, df_peaks, on='SampleName', how='inner')

# CSV file path
csv_file_frags = "frags_mqc.csv"
csv_file_peaks = "peaks_mqc.csv"

# Write header and DataFrame to the CSV file
with open(csv_file_frags, 'w') as file:
    # Write the header
    file.write(head_file_content_frags + "\n")
    # Write the DataFrame content
    df_frags.to_csv(file, index=False)

# Write header and DataFrame to the CSV file
with open(csv_file_peaks, 'w') as file:
    # Write the header
    file.write(head_file_content_peaks + "\n")
    # Write the DataFrame content
    df_peaks.to_csv(file, index=False)

print(f"CSV files successfully created at.")