#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q short.qb
#$ -t 1-22 -tc 16
#$ -N gnomAD_VEP_output_format
#$ -o /well/lindgren/George/Workflows/SS_constraint/Log/
#$ -e /well/lindgren/George/Workflows/SS_constraint/Log/

### Script that formats Ensembl VEP vcf output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:
# 1.CHROM, 2.POS, 3.ID, 4.REF, 5.ALT, 6.QUAL, 7.FILTER, 8.AF, 9.Gene, 10.Feature, 11.Feature_type, 12.Consequence, 
# 13.IMPACT, 14.SYMBOL, 15.SYMBOL_SOURCE, 16.BIOTYPE, 17.CANONICAL, 18.CCDS, 19.AF_male; 20.AF_female; 
# 21.non_neuro_AF_male; 22.non_neuro_AF_female; 23.controls_AF_male; 24.controls_AF_female; 25.InbreedingCoeff; 
# 26.VQSLOD; 27.n_alt_alleles


### $PED_ROOT is the rootname of your files 
PED_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/
#SGE_TASK_ID=22

#head -10000 gnomAD_VEP_output_chr17.vcf > gnomAD_VEP_output_chrdummy.vcf
	
### Save # lines to meta file
awk '/^#/' "$PED_ROOT"gnomAD_VEP_output_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"gnomAD_VEP_output_chr"$SGE_TASK_ID".meta

### Remove # lines
sed '/^#/ d' < "$PED_ROOT"gnomAD_VEP_output_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt

### Replace '|' with '\t to split INFO column
tr '|' '\t' < "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp"$SGE_TASK_ID"
mv "$PED_ROOT"tmp"$SGE_TASK_ID" "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt

### Replace ';' with '\t to split AC from INFO
tr ';' '\t' < "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp"$SGE_TASK_ID"
mv "$PED_ROOT"tmp"$SGE_TASK_ID" "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt

### remove "CSQ="
awk '{gsub("CSQ=", "");print}' "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp"$SGE_TASK_ID"
mv "$PED_ROOT"tmp"$SGE_TASK_ID" "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt

### Subset canonical protein-coding annotations
awk '$16 == "protein_coding" && $17 == "YES"' "$PED_ROOT"gnomAD_QCed_VEP_all_chr"$SGE_TASK_ID".txt > "$PED_ROOT"gnomAD_QCed_VEP_canPC_chr"$SGE_TASK_ID".txt
	


## summarise consequences
#awk -F '\t' '{print $12}' H_1000GP_QCed_VEP_canPC_dummy.txt | sort | uniq -c | sort -nr

## Count number of unique:
#cut -f 8 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Gene
#cut -f 9 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Feature
#cut -f 13 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # HGNC name
#cut -f 17 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # CCDD.ID 

### Check number of columns
#awk '{print NF}' H_1000GP_dummy_VEP.txt | sort -nu | tail -n 1





