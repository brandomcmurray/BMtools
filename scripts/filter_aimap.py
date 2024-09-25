#!/usr/bin/env python3

import os
import argparse

def parse_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith("Accession"):
                fields = line.strip().split('\t')
                data.append(fields)
    return data

def main():
    parser = argparse.ArgumentParser(
        description="Process and filter data from AIMAP output files.",
        epilog="© 2024 Brandon McMurray. All rights reserved."
    )
    parser.add_argument("-i", "--inputdir", default=".", help="specify the input file directory path")
    parser.add_argument("-x", "--infile_suffix", default="_full_result.txt", help="specify the suffix common to all input files")
    parser.add_argument("-o", "--outputdir", default=".", help="specify the path to save the output file")
    parser.add_argument("-n", "--outfile_name", default="output.txt", help="specify the name of the output file")
    parser.add_argument("-c", "--coverage", type=float, default=20, help="specify the Coverage cutoff (any numeric value greater than 1)")
    parser.add_argument("-e", "--edit_level", type=float, default=0.2, help="specify the Edit_level cutoff (any numeric value between 0 and 1)")
    parser.add_argument("--sort", action="store_true", help="sort output table based on Position and Sample fields")
    args = parser.parse_args()

    all_data = {}

    # Parsing all files and storing data in a dictionary
    for root, dirs, files in os.walk(os.path.expanduser(args.inputdir)):
        for file in files:
            if file.endswith(args.infile_suffix):
                sample_name = file.split("_")[0]
                file_path = os.path.join(root, file)
                all_data[sample_name] = parse_file(file_path)

    # Filtering and writing output
    with open(os.path.join(args.outputdir, args.outfile_name), 'w') as out_file:
        out_file.write('\t'.join(["Accession", "Position", "Old_base", "New_base", "Raw_read_depth",
                                  "Coverage", "Edit_level", "snp_coverage", "snp_f_coverage", "snp_r_coverage",
                                  "Gene_biotype", "Gene_name", "Gene_strand", "Product", "Amino_acid_change", "Sample"]) + '\n')

        output_data = []

        for sample_name, data in all_data.items():
            for entry in data:
                accession, position, old_base, new_base, raw_read_depth, coverage, edit_level = entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6]
                if float(coverage) >= args.coverage and float(edit_level) >= args.edit_level:
                    matching_entries = []
                    for other_sample, other_data in all_data.items():
                        if other_sample != sample_name:
                            for other_entry in other_data:
                                if entry[1:4] == other_entry[1:4]:
                                    matching_entries.append(other_sample)

                    if matching_entries:
                        output_data.append(entry + [sample_name])

        if args.sort:
            # Sort the output based on the "Position" and then "Sample" column
            sorted_output_data = sorted(output_data, key=lambda x: (int(x[1]), x[-1]))
        else:
            sorted_output_data = output_data

        # Write the sorted output to the file
        for entry in sorted_output_data:
            out_file.write('\t'.join(entry) + '\n')

if __name__ == "__main__":
    main()
