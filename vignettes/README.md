# Example analysis walkthroughs
## A-to-I RNA editing analysis walkthrough
### Filter AIMAP results for likely edits (-c/--coverage >20 reads and -e/--edit_level >20% alternative allele)
```
filter_aimap.py -c 20 -e 0.2 --sort --inputdir ~/aimap --outputdir ~/aimap -n filtered_aimap_full_result.txt
```

#### Set strand for every detection; A_G == "+" and T_C == "-"
```
awk -F'\t' 'BEGIN {OFS=FS} $11 == "CDS" && $3 == "T" && $4 == "C" { $13 = "-" } $11 == "CDS" && $3 == "A" && $4 == "G" { $13 = "+" } { print }' ~/aimap/filtered_aimap_full_result.txt > ~/aimap/filtered_aimap_full_result_updated.txt
```

#### Wrangle the data of the filtered AIMAP table into a non-redundant pseudo-BED file; include mean of edit_level, mean of snp_coverage, std.dev of snp_coverage, median of snp_coverage, and the file IDs
```
awk 'BEGIN {FS=OFS="\t"}
    NR > 1 {
        chrom = $1;
        start = $2;
        end = $2;
        name = $2"_"$3"_"$4;
        strand = ($11 == "CDS" ? $13 : "*");
        biotype = $11;
        gene = ($11 == "CDS" ? $12 : "NA");
        product = ($11 == "CDS" ? $14 : "NA");
        amino = ($11 == "CDS" ? $15 : "NA");
        value7 = $7;
        value8 = $8;
        lastField = $(NF);

        names[name]++;
        sum7[name] += value7;
        sum8[name] += value8;
        values8[name][names[name]] = value8;
        lastFields[name] = (lastFields[name] ? lastFields[name] "," : "") lastField;

        data[name] = chrom OFS start OFS end OFS name OFS "." OFS strand OFS biotype OFS gene OFS product OFS amino;
    }

    function median(arr, n) {
        if (n % 2) {
            return arr[int(n/2)+1];
        } else {
            return (arr[n/2] + arr[n/2+1]) / 2.0;
        }
    }

    END {
        for (name in names) {
            n = names[name];
            avg7 = sum7[name] / n;
            avg8 = sum8[name] / n;

            sumSquares = 0;
            for (i = 1; i <= n; i++) {
                sumSquares += (values8[name][i] - avg8)^2;
            }
            stddev8 = sqrt(sumSquares / n);

            asort(values8[name]);
            med8 = median(values8[name], n);

            print data[name], avg7, avg8, stddev8, med8, lastFields[name];
        }
    }' ~/aimap/filtered_aimap_full_result_updated.txt | sort -k2,2n > ~/aimap/aimap_nr_edit_data_full.txt
```

#### Subset the table to confident A-to-I modification detections in at least two biological replicates
```
awk '($4 ~ "A_G" || $4 ~ "T_C") && $NF ~ ","' ~/aimap/aimap_nr_edit_data_full.txt > ~/aimap/aimap_nr_edit_data.txt
```

#### Convert to BED; force the A-to-I modifications in RNA from intergenic regions to be stranded ("+" for A-to-G changes and "-" for T-to-C changes)
```
awk 'BEGIN {FS=OFS="\t"} { if ($4 ~ /A_G/) $6 = "+"; else if ($4 ~ /T_C/) $6 = "-"; print $1, $2, $3, $4, $5, $6 }' ~/aimap/aimap_nr_edit_data.txt > ~/aimap/aimap_nr_edit_data_stranded.bed
```

#### Convert to BED; subset the table to A-to-I modifications in mRNAs
```
awk 'BEGIN {FS=OFS="\t"} $7 == "CDS" {print $1, $2, $3, $4, $5, $6}' ~/aimap/aimap_nr_edit_data.txt > ~/aimap/aimap_nr_edit_data_CDS.bed
```

### Search for a common nucleotide pattern amongst editing sites
#### Find all genomic sites matching a specific nucleotide pattern
```
fuzznuc -sequence GCF_001886595.1_ASM188659v1_genomic.fna -pattern NAHG -pname TadA_motif -outfile ~/aimap/NAHG_genomic_positions_seqtable.fuzznuc -rformat seqtable -complement Y
fuzznuc -sequence GCF_001886595.1_ASM188659v1_genomic.fna -pattern NAHG -pname TadA_motif -outfile ~/aimap/NAHG_genomic_positions_excel.fuzznuc -rformat excel -complement Y
```

#### Merge the results of the two fuzznuc runs
```
seq_info_transfer.py -s ~/aimap/NAHG_genomic_positions_seqtable.fuzznuc -t ~/aimap/NAHG_genomic_positions_excel.fuzznuc -o ~/aimap/NAHG_genomic_positions_excel_updated.fuzznuc
```

#### Subset the updated.fuzznuc file to only those lines spanning positions specified in the aimap.bed file INCLUDING STRAND INFO
```
awk 'NR==FNR { ids[$1,$2,$3,$6]=$4; next } FNR==1 { print $0"\tID"; next } { for (id in ids) { split(id,a,SUBSEP); if ($1==a[1] && $2<=a[2] && $3>=a[3] && $4==a[4]) { print $0"\t"ids[id] } } }' ~/aimap/aimap_nr_edit_data_stranded.bed ~/aimap/NAHG_genomic_positions_excel_updated.fuzznuc > ~/aimap/NAHG_genomic_positions_excel_updated_subset_aimap_nr_edit_data_stranded.fuzznuc
```

#### Count only those sites corresponding to the edited "A" in the second position of the motif (this is a strand-aware method)
```
awk 'BEGIN {FS=OFS="\t"}
{
  split($8, id, "_")
  if (($4 == "+" && id[1] == $2 + 1) || ($4 == "-" && id[1] == $3 - 1)) {
    print
  }
}' ~/aimap/NAHG_genomic_positions_excel_updated_subset_aimap_nr_edit_data_stranded.fuzznuc | wc -l
```

### Compute local sequence secondary structures for all mRNAs with edits in them
#### Extract 117 nts (-58,0,+58) centred on edit site and find negative control sequences
```
findMotifsGenome.pl ~/aimap/aimap_nr_edit_data_stranded.bed GCF_001886595.1_ASM188659v1_genomic.fna ~/aimap/homer117_all_rna_fasta -size -58,58 -len 6 -rna -dumpFasta
```
```
findMotifsGenome.pl ~/aimap/aimap_nr_edit_data_CDS.bed GCF_001886595.1_ASM188659v1_genomic.fna ~/aimap/homer117_CDS_rna_fasta -size -58,58 -len 6 -rna -dumpFasta
```

#### Extract 89 nts for exactly the 73-nts of the positive control tRNA-Arg[ACG] (N.B. edit site at position 4090444; -33,0,+39)
```
samtools faidx GCF_001886595.1_ASM188659v1_genomic.fna NZ_CP018074.1:4090403-4090491 > ~/aimap/trna_arg.fa
```

#### Convert dumped fasta target and background sequences from DNA to RNA
For all edit positions
```
awk '/^[^>]/{gsub("T","U")}1' ~/aimap/homer117_all_rna_fasta/target.fa > ~/aimap/RNAfold/target_117_all_rna.fa
awk '/^[^>]/{gsub("T","U")}1' ~/aimap/homer117_all_rna_fasta/background.fa > ~/aimap/RNAfold/background_117_all_rna.fa
```
For only mRNA edit positions
```
awk '/^[^>]/{gsub("T","U")}1' ~/aimap/homer117_CDS_rna_fasta/target.fa > ~/aimap/RNAfold/target_117_cds_rna.fa
awk '/^[^>]/{gsub("T","U")}1' ~/aimap/homer117_CDS_rna_fasta/background.fa > ~/aimap/RNAfold/background_117_cds_rna.fa
```
For the positive tRNA-Arg transcript
```
awk '/^[^>]/{gsub("T","U")}1' ~/aimap/trna_arg.fa > ~/aimap/RNAfold/positive_trna_arg.fa
```

#### Calculate the minimum free energy of every 17-nt sliding window (length of anticodon arm of tRNA) in a 117-nt sequence centred on an edit site
```
sliding_window_mfe_calculator.sh target_117_all_rna.fa ~/aimap/RNAfold/100nt_transcripts/target_100pos_all_rna_mfe_results.txt 17
sliding_window_mfe_calculator.sh background_117_all_rna.fa ~/aimap/RNAfold/100nt_transcripts/background_100pos_all_rna_mfe_results.txt 17
sliding_window_mfe_calculator.sh target_117_cds_rna.fa ~/aimap/RNAfold/100nt_transcripts/target_100pos_cds_rna_mfe_results.txt 17
sliding_window_mfe_calculator.sh background_117_cds_rna.fa ~/aimap/RNAfold/100nt_transcripts/background_100pos_cds_rna_mfe_results.txt 17
```
Run the program for the positive control, but then also update the Index/Position column for offset edit position
```
sliding_window_mfe_calculator.sh ~/aimap/RNAfold/positive_trna_arg.fa ~/aimap/RNAfold/100nt_transcripts/positive_trna_arg_mfe_results_temp.txt 17
```
```
awk 'BEGIN {OFS="\t"} NR==1 {print; next} {$4 = -33 + (NR-2); print}' ~/aimap/RNAfold/100nt_transcripts/positive_trna_arg_mfe_results_temp.txt > ~/aimap/RNAfold/100nt_transcripts/positive_trna_arg_mfe_results.txt; rm -f ~/aimap/RNAfold/100nt_transcripts/positive_trna_arg_mfe_results_temp.txt
```

#### Calculate the mean minimum free energy and standard deviation of each nucleotide position
```
for file in target_100pos_all_rna background_100pos_all_rna target_100pos_cds_rna background_100pos_cds_rna; do
awk 'BEGIN {FS=OFS="\t"} 
    NR > 1 {
    idx[$4] = $4;
    mfe[$4] += $3;
    count[$4]++;
    mfe_squared[$4] += ($3)^2;
} END {
    print "Index\tmean_mfe\tsd_mfe";
    for (i in idx) {
        mean = mfe[i] / count[i];
        variance = (mfe_squared[i] / count[i]) - (mean^2);
        sd = sqrt(variance);
        printf "%s\t%.5f\t%.5f\n", i, mean, sd;
    }
}' ~/aimap/RNAfold/100nt_transcripts/$file\_mfe_results.txt | { read header; echo "$header"; sort -k1,1n; } > ~/aimap/RNAfold/100nt_transcripts/$file\_mfe_results_mean_sd.txt
done
```
