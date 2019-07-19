#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -N ClinVar_VEP_output_format
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that formats Ensembl VEP output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:


### set rootnames of your files
IN_ROOT=/well/lindgren/George/Data/ClinVar/formatted/
OUT_ROOT=/well/lindgren/George/Data/ClinVar/formatted/

### set file names
IN_FILE=ClinVar_pathogenic_snps_QCed_VEPout.vcf
OUT_FILE=ClinVar_pathogenic_snps_QCed_VEP_v94.vcf

### Save # lines 
awk '/^#/' "$IN_ROOT""$IN_FILE" > "$IN_ROOT""$IN_FILE".meta

### Remove # lines
sed '/^#/ d' < "$IN_ROOT""$IN_FILE" > "$IN_ROOT"tmp_file

### Replace '|' with '\t to split INFO column
tr '|' '\t' < "$IN_ROOT"tmp_file > "$IN_ROOT"tmp
mv "$IN_ROOT"tmp "$IN_ROOT"tmp_file

### remove "CSQ="
awk '{gsub("CSQ=", "");print}' "$IN_ROOT"tmp_file > "$IN_ROOT"tmp && mv "$IN_ROOT"tmp "$IN_ROOT"tmp_file

### Subset canonical protein-coding annotations
#awk '$15 == "protein_coding" && $16 == "YES"' ClinVar_VEP_output.txt > ClinVar_VEP_canPC.txt

### rename file
mv "$IN_ROOT"tmp_file "$OUT_ROOT""$OUT_FILE"