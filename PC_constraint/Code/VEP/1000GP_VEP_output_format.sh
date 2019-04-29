#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -t 1-22 -tc 16
#$ -N thouGP_VEP_output_format
#$ -o /well/lindgren/George/Workflows/PC_constraint/Log/
#$ -e /well/lindgren/George/Workflows/PC_constraint/Log/

### Script that formats Ensembl VEP vcf output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:
# 1.CHROM, 2.POS, 3.ID, 4.REF, 5.ALT, 6.QUAL, 7.FILTER, 8.AC, 9.Gene, 10.Feature, 11.Feature_type, 12.Consequence, 
# 13.IMPACT, 14.SYMBOL, 15.SYMBOL_SOURCE, 16.BIOTYPE, 17.CANONICAL, 18.CCDS, 19.AF


### $PED_ROOT is the rootname of your files 
PED_ROOT=/well/lindgren/George/Data/1000GP/vcf_QCed_VEP/
# SGE_TASK_ID=dummy
# head -10000 1000GP_all_VEP_v94_output_chr1.vcf > 1000GP_all_VEP_v94_output_chr$SGE_TASK_ID.vcf

	echo "Processing 1000GP_all_VEP_v94_output_chr$SGE_TASK_ID.vcf"
	
	### Save # lines to meta file
	awk '/^#/' "$PED_ROOT"1000GP_all_VEP_v94_output_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"1000GP_all_VEP_v94_output_chr"$SGE_TASK_ID".meta

	### Remove # lines
	sed '/^#/ d' < "$PED_ROOT"1000GP_all_VEP_v94_output_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt

	### Replace '|' with '\t to split INFO column
	tr '|' '\t' < "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt

	### Replace ';' with '\t to split AC from INFO
	tr ';' '\t' < "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt

	### remove "CSQ="
	awk '{gsub("CSQ=", "");print}' "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt

	### Subset canonical protein-coding annotations
	awk '$16 == "protein_coding" && $17 == "YES"' "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt >> "$PED_ROOT"1000GP_phase3_QCed_VEP_v94_canPC_PASS.txt
	
	echo "1000GP_all_VEP_v94_output_chr$SGE_TASK_ID.vcf complete!"


## summarise consequences
#awk -F '\t' '{print $12}' H_1000GP_QCed_VEP_canPC_dummy.txt | sort | uniq -c | sort -nr

## Count number of unique:
#cut -f 8 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Gene
#cut -f 9 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Feature
#cut -f 13 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # HGNC name
#cut -f 17 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # CCDD.ID 

### Check number of columns
#awk '{print NF}' H_1000GP_dummy_VEP.txt | sort -nu | tail -n 1





