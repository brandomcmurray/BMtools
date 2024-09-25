#!/usr/bin/env python3

import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Extract sequences from a multi-FASTA file based on a list of locus_tags.')
    parser.add_argument('-l', '--list', required=True, help='Path to the input locus_tag list file')
    parser.add_argument('-f', '--fasta', required=True, help='Path to the multi-FASTA file')
    parser.add_argument('-o', '--out', required=True, help='Path to the output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print unmatched locus_tags to the terminal')
    return parser.parse_args()

def extract_sequences(locus_tags, fasta_file):
    sequences = {}
    matched_tags = set()
    with open(fasta_file, 'r') as f:
        header = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence:
                    locus_tag = header.split('locus_tag=')[1].split()[0]
                    if locus_tag in locus_tags:
                        sequences[header] = ''.join(sequence)
                        matched_tags.add(locus_tag)
                header = line
                sequence = []
            else:
                sequence.append(line)

        if header and sequence:
            locus_tag = header.split('locus_tag=')[1].split()[0]
            if locus_tag in locus_tags:
                sequences[header] = ''.join(sequence)
                matched_tags.add(locus_tag)

    return sequences, matched_tags

def write_output(sequences, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as out:
        for header, seq in sequences.items():
            out.write(f"{header}\n{seq}\n")

def main():
    args = parse_args()

    # Read the locus tags from the input file
    with open(args.list, 'r') as f:
        locus_tags = {line.strip() for line in f if line.strip()}

    # Extract sequences from the FASTA file
    sequences, matched_tags = extract_sequences(locus_tags, args.fasta)

    # Write the extracted sequences to the output file
    write_output(sequences, args.out)

    # Print unmatched locus_tags if verbose is enabled
    if args.verbose:
        unmatched_tags = locus_tags - matched_tags
        if unmatched_tags:
            print("Unmatched locus_tags:")
            for tag in unmatched_tags:
                print(tag)
        else:
            print("All locus_tags matched.")

if __name__ == '__main__':
    main()
