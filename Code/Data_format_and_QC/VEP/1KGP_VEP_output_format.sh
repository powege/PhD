#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -t 1-22 -tc 16
#$ -N thouGP_VEP_output_format
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that formats Ensembl VEP vcf output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:
# 1.CHROM, 2.POS, 3.ID, 4.REF, 5.ALT, 6.QUAL, 7.FILTER, 8.AC, 9.Gene, 10.Feature, 11.Feature_type, 12.Consequence, 
# 13.IMPACT, 14.SYMBOL, 15.SYMBOL_SOURCE, 16.BIOTYPE, 17.CANONICAL, 18.CCDS, 19.AF


### set file root
PED_ROOT=/well/lindgren/George/Data/1KGP/Variants/vcf_QCed_VEP/

### set file name
IN_FILE=1KGP_phase3_snps_QCed_VEPout_all_chr
OUT_FILE=1KGP_phase3_snps_QCed_VEP_v94_chr
CAN_PC_FILE=1KGP_phase3_snps_QCed_VEP_v94_canPC_chr

# SGE_TASK_ID=dummy
# head -10000 "$PED_ROOT"1KGP_phase3_snps_QCed_VEPout_all_chr1.vcf > "$PED_ROOT"1KGP_phase3_snps_QCed_VEPout_all_chr"$SGE_TASK_ID".vcf

	echo "Processing $IN_FILE$SGE_TASK_ID.vcf"
	
	### Save # lines to meta file
	awk '/^#/' "$PED_ROOT""$IN_FILE""$SGE_TASK_ID".vcf > "$PED_ROOT""$IN_FILE""$SGE_TASK_ID".meta

	### Remove # lines
	sed '/^#/ d' < "$PED_ROOT""$IN_FILE""$SGE_TASK_ID".vcf > "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt

	### Replace '|' with '\t' to split INFO column
	tr '|' '\t' < "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt

	### Replace ';' with '\t' to split AC from INFO
	tr ';' '\t' < "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt

	### remove "CSQ="
	awk '{gsub("CSQ=", "");print}' "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt

	### Replace '\t' with ','
	# tr '\t' ',' < "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	# mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt

	### Subset canonical protein-coding annotations
	awk '$16 == "protein_coding" && $17 == "YES"' "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt > "$PED_ROOT""$CAN_PC_FILE""$SGE_TASK_ID".vcf
	
	### rename output
	mv "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt "$PED_ROOT""$OUT_FILE""$SGE_TASK_ID".vcf
	
	### remove temp files
	# rm "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt
	
	echo "$IN_FILE$SGE_TASK_ID.vcf complete!"


## summarise consequences
#awk -F '\t' '{print $12}' H_1000GP_QCed_VEP_canPC_dummy.txt | sort | uniq -c | sort -nr

## Count number of unique:
#cut -f 8 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Gene
#cut -f 9 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # Feature
#cut -f 13 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # HGNC name
#cut -f 17 H_1000GP_QCed_canPC_dummy.txt | sort | uniq | wc -l # CCDD.ID 

### Check number of columns
#awk '{print NF}' "$PED_ROOT"1KGP_phase3_snps_QCed_VEP_v94_chr"$SGE_TASK_ID".txt | sort -nu | tail -n 1





