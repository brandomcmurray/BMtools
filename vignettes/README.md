# Example analysis walkthroughs

## A-to-I RNA editing analysis walkthrough

### Filtering AIMAP results for likely edits
```
./filter_aimap.py -c 20 -e 0.2 --sort --inputdir ~/results/rnaseq/aimap --outputdir ~/results/rnaseq/aimap -n filtered_aimap_full_result.txt
```

#### Set strand for every detection; A_G == "+" and T_C == "-"
```
awk -F'\t' 'BEGIN {OFS=FS} $11 == "CDS" && $3 == "T" && $4 == "C" { $13 = "-" } $11 == "CDS" && $3 == "A" && $4 == "G" { $13 = "+" } { print }' ~/results/rnaseq/aimap/filtered_aimap_full_result.txt > ~/results/rnaseq/aimap/filtered_aimap_full_result_updated.txt
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

            # Calculate standard deviation for field 8
            sumSquares = 0;
            for (i = 1; i <= n; i++) {
                sumSquares += (values8[name][i] - avg8)^2;
            }
            stddev8 = sqrt(sumSquares / n);

            # Calculate standard error for field 8
            #stderror8 = stddev8 / sqrt(n);

            # Calculate median for field 8
            asort(values8[name]);
            med8 = median(values8[name], n);

            print data[name], avg7, avg8, stddev8, med8, lastFields[name];
        }
    }' ~/results/rnaseq/aimap/filtered_aimap_full_result_updated.txt | sort -k2,2n > ~/results/rnaseq/aimap/aimap_nr_edit_data_full.txt
```

#### Subset the table to confident A-to-I modification detections in at least two biological replicates
```
awk '($4 ~ "A_G" || $4 ~ "T_C") && $NF ~ ","' ~/results/rnaseq/aimap/aimap_nr_edit_data_full.txt > ~/results/rnaseq/aimap/aimap_nr_edit_data.txt
```

#### Convert to BED; force the A-to-I modifications in RNA from intergenic regions to be stranded ("+" for A-to-G changes and "-" for T-to-C changes)
```
awk 'BEGIN {FS=OFS="\t"} { if ($4 ~ /A_G/) $6 = "+"; else if ($4 ~ /T_C/) $6 = "-"; print $1, $2, $3, $4, $5, $6 }' ~/results/rnaseq/aimap/aimap_nr_edit_data.txt > ~/results/rnaseq/aimap/aimap_nr_edit_data_stranded.bed
```

#### Convert to BED; subset the table to A-to-I modifications in mRNAs
```
awk 'BEGIN {FS=OFS="\t"} $7 == "CDS" {print $1, $2, $3, $4, $5, $6}' ~/results/rnaseq/aimap/aimap_nr_edit_data.txt > ~/results/rnaseq/aimap/aimap_nr_edit_data_CDS.bed
```

### Search for a common nucleotide pattern amongst editing sites
#### Find all genomic sites matching a specific nucleotide pattern
```
fuzznuc -sequence ~/test/genome/GCF_001886595.1_ASM188659v1_genomic.fna -pattern NAHG -pname TadA_motif -outfile ~/test/NAHG_genomic_positions_seqtable.fuzznuc -rformat seqtable -complement Y
fuzznuc -sequence ~/test/genome/GCF_001886595.1_ASM188659v1_genomic.fna -pattern NAHG -pname TadA_motif -outfile ~/test/NAHG_genomic_positions_excel.fuzznuc -rformat excel -complement Y
```

#### Merge the results of the two fuzznuc runs
```
seq_info_transfer.py -s ~/test/NAHG_genomic_positions_seqtable.fuzznuc -t ~/test/NAHG_genomic_positions_excel.fuzznuc -o ~/test/NAHG_genomic_positions_excel_updated.fuzznuc
```

#### Subset the updated.fuzznuc file to only those lines spanning positions specified in the aimap.bed file INCLUDING STRAND INFO
```
awk 'NR==FNR { ids[$1,$2,$3,$6]=$4; next } FNR==1 { print $0"\tID"; next } { for (id in ids) { split(id,a,SUBSEP); if ($1==a[1] && $2<=a[2] && $3>=a[3] && $4==a[4]) { print $0"\t"ids[id] } } }' ~/results/rnaseq/aimap/aimap_nr_edit_data_stranded.bed ~/test/NAHG_genomic_positions_excel_updated.fuzznuc > ~/results/rnaseq/aimap/NAHG_genomic_positions_excel_updated_subset_aimap_nr_edit_data_stranded.fuzznuc
```

#### Count only those sites corresponding to the edited "A" in the second position of the motif (this is a strand-aware method)
awk 'BEGIN {FS=OFS="\t"}
{
  split($8, id, "_")
  if (($4 == "+" && id[1] == $2 + 1) || ($4 == "-" && id[1] == $3 - 1)) {
    print
  }
}' ~/results/rnaseq/aimap/NAHG_genomic_positions_excel_updated_subset_aimap_nr_edit_data_stranded.fuzznuc | wc -l
