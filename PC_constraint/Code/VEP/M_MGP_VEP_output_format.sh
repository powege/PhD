#!/bin/bash

### Script that formats Ensembl VEP output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:
# 1.CHROM, 2.POS, 3.ID, 4.REF, 5.ALT, 6.QUAL, 7.FILTER, 8.Gene, 9.Feature, 10.Feature_type, 11.Consequence, 
# 12.IMPACT, 13.SYMBOL, 14.SYMBOL_SOURCE, 15.BIOTYPE, 16.CANONICAL, 17.CCDS,

### $PED_ROOT is the rootname of your files 
#PED_ROOT=/Users/g.powell/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/
cd /Users/g.powell/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP

### Save # lines 
awk '/^#/' M_MGP_QCed_VEP_output_all.vcf > M_MGP_QCed_VEP_output_all.meta

### Remove # lines
sed '/^#/ d' < M_MGP_QCed_VEP_output_all.vcf > M_MGP_QCed_VEP_all.txt

### Replace '|' with '\t to split INFO column
tr '|' '\t' < M_MGP_QCed_VEP_all.txt > tmp
mv tmp M_MGP_QCed_VEP_all.txt

### remove "CSQ="
awk '{gsub("CSQ=", "");print}' M_MGP_QCed_VEP_all.txt > tmp && mv tmp M_MGP_QCed_VEP_all.txt

### Subset canonical protein-coding annotations
awk '$15 == "protein_coding" && $16 == "YES"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_canPC_all.txt

### split by chromosome
awk '$1 == "1"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr1.txt
awk '$1 == "2"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr2.txt
awk '$1 == "3"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr3.txt
awk '$1 == "4"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr4.txt
awk '$1 == "5"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr5.txt
awk '$1 == "6"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr6.txt
awk '$1 == "7"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr7.txt
awk '$1 == "8"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr8.txt
awk '$1 == "9"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr9.txt
awk '$1 == "10"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr10.txt
awk '$1 == "11"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr11.txt
awk '$1 == "12"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr12.txt
awk '$1 == "13"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr13.txt
awk '$1 == "14"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr14.txt
awk '$1 == "15"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr15.txt
awk '$1 == "16"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr16.txt
awk '$1 == "17"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr17.txt
awk '$1 == "18"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr18.txt
awk '$1 == "19"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chr19.txt
awk '$1 == "X"' M_MGP_QCed_VEP_all.txt > M_MGP_QCed_VEP_all_chrX.txt


## summarise consequences
#awk -F '\t' '{print $11}' M_MGP_QCed_VEP_canPC_all.txt | sort | uniq -c | sort -nr

## Count number of unique:
#cut -f 8 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # Gene
#cut -f 9 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # Feature
#cut -f 13 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # HGNC name
#cut -f 17 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # CCDD.ID 

### Check number of columns
#awk '{print NF}' M_MGP_QCed_VEP_canPC_all.txt | sort -nu | tail -n 1

