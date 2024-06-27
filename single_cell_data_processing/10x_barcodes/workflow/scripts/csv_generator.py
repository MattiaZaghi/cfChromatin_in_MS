import csv
import argparse

def main():
    # Get input arguments from command line
    args = parse_arguments()

    # Create a list to store CSV rows
    csv_rows = []

    # Loop through the input directories and create rows for the CSV
    sample = args.sample_data
    fragments = args.fragments_path
    single_cells=args.cells_path
    csv_rows.append([sample,fragments, single_cells])

    # Write data to CSV file
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['library_id', 'fragments', 'cells'])
        csvwriter.writerows(csv_rows)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate CSV file for epigenomic modality and RNA samples')
    parser.add_argument('--sample_data', help='Sample data', required=True)
    parser.add_argument('--fragments_path', help='fragments tsv file', required=True)
    parser.add_argument('--cells_path', help='single cells csv file', required=True)
    parser.add_argument('-o', '--output', help='output CSV file', required=True)

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()
