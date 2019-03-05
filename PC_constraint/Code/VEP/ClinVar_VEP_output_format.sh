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
awk '/^#/' ClinVar_VEP_output.vcf > ClinVar_VEP_output.meta

### Remove # lines
sed '/^#/ d' < ClinVar_VEP_output.vcf > ClinVar_VEP_output.txt

### Replace '|' with '\t to split INFO column
tr '|' '\t' < ClinVar_VEP_output.txt > tmp
mv tmp ClinVar_VEP_output.txt

### remove "CSQ="
awk '{gsub("CSQ=", "");print}' ClinVar_VEP_output.txt > tmp && mv tmp ClinVar_VEP_output.txt

### Subset canonical protein-coding annotations
#awk '$15 == "protein_coding" && $16 == "YES"' ClinVar_VEP_output.txt > ClinVar_VEP_canPC.txt
