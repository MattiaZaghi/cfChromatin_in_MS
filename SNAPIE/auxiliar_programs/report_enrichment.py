
#!/usr/bin/env python
# coding: utf-8
import os
import pandas as pd
import argparse

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
        print(f"Error: The file '{{str_path}}' was not found.")
        return ""
    except IOError as e:
        print(f"Error: Could not read the file '{{str_path}}': {e}")
        return ""
    return content


def load_files_to_dataframe(file_list, mark,SampleName):
    """
    Loads the content of multiple space-separated files into a DataFrame.
    Assumes the first line of each file contains column names.

    :param file_list: List of file paths to read.
    :param mark: The mark parameter to filter the DataFrame.
    :return: A concatenated pandas DataFrame containing the data from all files.
    """
    dataframes = []

    for file_path in file_list:
        try:
            # Load the file into a DataFrame
            df = pd.read_csv(file_path, delim_whitespace=True, header=0)
            df['File Name'] = file_path  # Add a column with the file name
            dataframes.append(df)
        except FileNotFoundError:
            print(f"Error: The file '{{file_path}}' was not found.")
        except Exception as e:
            print(f"Error: Could not process the file '{{file_path}}': {e}")

    # Concatenate all DataFrames
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df['mark_upper'] = combined_df['mark'].str.upper()
        combined_df['name_upper'] = combined_df['name'].str.upper()
        combined_df = combined_df[(combined_df['mark_upper'] == mark.upper()) & (combined_df['name_upper'] == SampleName.upper())]
        combined_df = combined_df.drop(columns=['mark_upper','name_upper', 'File Name'], errors='ignore')
        combined_df = combined_df.rename(columns={'name': 'SampleName'})

        return combined_df
    else:
        print("No valid files were provided.")
        return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(description="Process enrichment files with a specified mark.")
    parser.add_argument("--mark", required=True, help="Specify the mark to filter the enrichment files.")
    parser.add_argument("--samplename", required=True, help="Specify the SampleName to filter the enrichment files.")
    args = parser.parse_args()

    mark = args.mark
    samplename = args.samplename

    current_directory = os.getcwd()

    # Find and process enrichment files
    suffix_enrichment = '_enrichment_states.csv'
    enrichment_files = find_files_with_suffix(current_directory, suffix_enrichment)
    df_enrichment = load_files_to_dataframe(enrichment_files, mark, samplename)

    # Check if the DataFrame is empty after filtering
    if df_enrichment.empty:
        print(f"No records found for the specified mark '{mark}'. File creation skipped.")
        return

    # CSV file path
    csv_file = f"enrichment_{mark.upper()}_{samplename.upper()}_report.csv"
        
    # Write the DataFrame content
    df_enrichment.to_csv(csv_file, index=False)

    print(f"CSV file successfully created at '{{csv_file}}'.")


if __name__ == "__main__":
    main()
