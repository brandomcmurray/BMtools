#!/usr/bin/env python3

import argparse
import re

def parse_fuzznuc(file_path):
    sequences = {}
    current_sequence = None

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("# Sequence:"):
                current_sequence = re.search(r'NZ_CP0180\d+\.\d+', line).group()
            elif line.strip() and not line.startswith("#"):
                fields = line.split()
                key = (current_sequence, fields[0], fields[1], fields[2], fields[3], fields[4])
                sequences[key] = fields[5]

    return sequences

def parse_txt(file_path):
    entries = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith("SeqName"):
                fields = line.split()
                key = (fields[0], fields[1], fields[2], fields[4], fields[5], fields[6])
                entries.append(key)

    return entries

def write_output(matches, output_file_path):
    with open(output_file_path, 'w') as output_file:
        header = "SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\tSequence\n"
        output_file.write(header)

        for entry, sequence in matches.items():
            line = '\t'.join(entry) + f'\t{sequence}\n'
            output_file.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Copy the nucleotide sequence data from the EMBOSS fuzznuc seqtable output to the excel output.",
        epilog="© 2024 Brandon McMurray. All rights reserved."
    )
    parser.add_argument("-s", "--seqtable", default="./motif_seqtable.fuzznuc", help="specify the path to the fuzznuc seqtable file")
    parser.add_argument("-t", "--table", default="./motif_excel.fuzznuc", help="specify the path to the fuzznuc excel table file")
    parser.add_argument("-o", "--output", default="./motif_excel_updated.fuzznuc", help="specify the path to save the updated fuzznuc excel table file")
    args = parser.parse_args()

    fuzznuc_sequences = parse_fuzznuc(args.seqtable)
    txt_entries = parse_txt(args.table)

    matches = {}

    for entry in txt_entries:
        if entry in fuzznuc_sequences:
            sequence = fuzznuc_sequences[entry]
            matches[entry] = sequence

    write_output(matches, args.output)
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
