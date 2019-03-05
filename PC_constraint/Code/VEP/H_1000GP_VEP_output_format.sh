#!/bin/bash

### Script that formats Ensembl VEP vcf output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:
# 1.CHROM, 2.POS, 3.ID, 4.REF, 5.ALT, 6.QUAL, 7.FILTER, 8.AC, 9.Gene, 10.Feature, 11.Feature_type, 12.Consequence, 
# 13.IMPACT, 14.SYMBOL, 15.SYMBOL_SOURCE, 16.BIOTYPE, 17.CANONICAL, 18.CCDS, 19.AF,
# 20.EAS_AF, 21.AMR_AF, 22.AFR_AF, 23.EUR_AF, 24.SAS_AF


### $PED_ROOT is the rootname of your files 
PED_ROOT=/Users/g.powell/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/
#cd /Users/g.powell/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do 
	echo "Processing H_1000GP_QCed_VEP_output_all_chr$CHR.vcf"
	
	### Save # lines to meta file
	awk '/^#/' "$PED_ROOT"H_1000GP_QCed_VEP_output_all_chr"$CHR".vcf > "$PED_ROOT"H_1000GP_QCed_VEP_output_all_chr"$CHR".meta

	### Remove # lines
	sed '/^#/ d' < "$PED_ROOT"H_1000GP_QCed_VEP_output_all_chr"$CHR".vcf > "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt

	### Replace '|' with '\t to split INFO column
	tr '|' '\t' < "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt > "$PED_ROOT"tmp
	mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt

	### Replace ';' with '\t to split AC from INFO
	tr ';' '\t' < "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt > "$PED_ROOT"tmp
	mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt

	### remove "CSQ="
	awk '{gsub("CSQ=", "");print}' "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt

	### Subset canonical protein-coding annotations
	awk '$16 == "protein_coding" && $17 == "YES"' "$PED_ROOT"H_1000GP_QCed_VEP_all_chr"$CHR".txt > "$PED_ROOT"H_1000GP_QCed_VEP_all_canPC_chr"$CHR".txt
	
	echo "H_1000GP_QCed_VEP_output_all_chr$CHR.vcf complete!"
done


## summarise consequences
#awk -F '\t' '{print $12}' H_1000GP_QCed_VEP_canPC_dummy.txt | sort | uniq -c | sort -nr

## Count number of unique:
#cut -f 8 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Gene
#cut -f 9 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Feature
#cut -f 13 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # HGNC name
#cut -f 17 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # CCDD.ID 

### Check number of columns
#awk '{print NF}' H_1000GP_dummy_VEP.txt | sort -nu | tail -n 1





