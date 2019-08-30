#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -t 1-19 -tc 16
#$ -N MGP_VEP_output_format
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that formats Ensembl VEP output
# outputs: tab delimited file with all variants; tab delimited file with protein coding annotations
# Columns:
# 1.CHROM, 2.POS, 3.ID, 4.REF, 5.ALT, 6.QUAL, 7.FILTER, 8.Gene, 9.Feature, 10.Feature_type, 11.Consequence, 
# 12.IMPACT, 13.SYMBOL, 14.SYMBOL_SOURCE, 15.BIOTYPE, 16.CANONICAL, 17.CCDS,

### set file root
PED_ROOT=/well/lindgren/George/Data/MGP/Variants/vcf_QCed_VEP/

### set file name
# IN_FILE=MGP_v5_allMUSMUS_snps_QCed_VEPout_hom_chr
# OUT_FILE=MGP_v5_allMUSMUS_snps_QCed_VEP_v94_chr
IN_FILE=MGP_v5_allSTRAIN_snps_QCed_VEPout_hom_chr
OUT_FILE=MGP_v5_allSTRAIN_snps_QCed_VEP_v94_chr
CAP_PC_FILE=MGP_v5_allSTRAIN_snps_QCed_VEP_v94_canPC_chr

# SGE_TASK_ID=dummy
# head -10000 MGP_v5_allMUSMUS_snps_QCed_VEPout_hom_chr1.vcf > MGP_v5_allMUSMUS_snps_QCed_VEPout_hom_chr$SGE_TASK_ID.vcf

	echo "Processing $IN_FILE$SGE_TASK_ID.vcf"
	
	### Save # lines to meta file
	awk '/^#/' "$PED_ROOT""$IN_FILE""$SGE_TASK_ID".vcf > "$PED_ROOT""$IN_FILE""$SGE_TASK_ID".meta

	### Remove # lines
	sed '/^#/ d' < "$PED_ROOT""$IN_FILE""$SGE_TASK_ID".vcf > "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt

	### Replace '|' with '\t' to split INFO column
	tr '|' '\t' < "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt

	### remove "CSQ="
	awk '{gsub("CSQ=", "");print}' "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt

  ### Replace '\t' with ','
#   tr '\t' ',' < "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_"$SGE_TASK_ID"
# 	mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt

	### Subset canonical protein-coding annotations
	awk '$15 == "protein_coding" && $16 == "YES"' "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt > "$PED_ROOT""$CAP_PC_FILE""$SGE_TASK_ID".vcf
	
	### rename output
	mv "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt "$PED_ROOT""$OUT_FILE""$SGE_TASK_ID".vcf
	
	### remove temp files
	# rm "$PED_ROOT"MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt
	
	echo "$IN_FILE$SGE_TASK_ID.vcf complete!"




## summarise consequences
#awk -F '\t' '{print $11}' M_MGP_QCed_VEP_canPC_all.txt | sort | uniq -c | sort -nr

## Count number of unique:
#cut -f 8 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # Gene
#cut -f 9 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # Feature
#cut -f 13 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # HGNC name
#cut -f 17 M_MGP_QCed_VEP_canPC_all.txt | sort | uniq | wc -l # CCDD.ID 

### Check number of columns
#awk '{print NF}' M_MGP_QCed_VEP_canPC_all.txt | sort -nu | tail -n 1

