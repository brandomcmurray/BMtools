#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_fasta> <output_file> <window_size>"
    exit 1
fi

input_fasta=$1
output_file=$2
window_size=$3

if ! [[ "$window_size" =~ ^[0-9]+$ ]] || [ "$window_size" -lt 1 ]; then
    echo "Error: Window size must be a positive integer."
    exit 1
fi

echo -e "Header\tSequence\tMFE\tIndex" > "$output_file"

header=""
sequence=""
while read -r line; do
    if [ "${line:0:1}" = ">" ]; then
        if [ -n "$header" ]; then
            seq_length=${#sequence}

            for ((i=0; i<seq_length-window_size+1; i++)); do
                window=${sequence:i:window_size}
                mfe=$(echo -e "$window\n" | /usr/local-centos6/ViennaRNA/Progs/RNAfold -T30 --noPS | awk 'NR==2{print $NF}' | tr -d '()')
                index=$((i - (seq_length - window_size) / 2))
                echo -e "${header:1}\t$window\t$mfe\t$index" >> "$output_file"
            done
        fi

        header="$line"
        sequence=""
    else
        sequence+="$line"
    fi
done < "$input_fasta"

if [ -n "$header" ] && [ -n "$sequence" ]; then
    seq_length=${#sequence}

    for ((i=0; i<seq_length-window_size+1; i++)); do
        window=${sequence:i:window_size}
        mfe=$(echo -e "$window\n" | /usr/local-centos6/ViennaRNA/Progs/RNAfold -T30 --noPS | awk 'NR==2{print $NF}' | tr -d '()')
        index=$((i - (seq_length - window_size) / 2))
        echo -e "${header:1}\t$window\t$mfe\t$index" >> "$output_file"
    done
fi
