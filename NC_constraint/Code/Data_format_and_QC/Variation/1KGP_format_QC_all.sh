#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -t 1-22 -tc 16
#$ -N thouGP_vcf_QC
#$ -o /well/lindgren/George/Workflows/NC_constraint/Log/
#$ -e /well/lindgren/George/Workflows/NC_constraint/Log/

### Script that QCs 1000GP vcfs by chromosome.
### QC: PASS filter status; SNV 
### Output cols: 1.CHROM 2.POS 3.ID 4.REF 5.ALT 6.QUAL 7.FILTER 8.AC 9.AF 

### Set working directory
PED_ROOT=/well/lindgren/George/Data/1000GP/vcf_raw/
OUT_ROOT=/well/lindgren/George/Data/1000GP/vcf_QCed_VEP/

# SGE_TASK_ID=dummy
# head -10000 "$PED_ROOT"ALL.chr6_GRCh38_sites.20170504.vcf > "$PED_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38_sites.20170504.vcf

	echo "Processing Chromosome $SGE_TASK_ID"

	# Save # lines 
	awk '/^#/' "$PED_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38_sites.20170504.vcf  > "$PED_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38_sites.20170504.vcf.meta
	# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER INFO
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' "$PED_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38_sites.20170504.vcf > "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	# Subset variants with "PASS" filter status (this removes header and meta lines)
	awk '$7 == "PASS"' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	# Subset variants with alternate base == 1
	awk 'length($5)<=1' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	# Subset variants with reference base == 1
	awk 'length($4)<=1' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	### subset variants with REF == A, T, G, or C
	awk '$4=="A" || $4=="T" || $4=="G" || $4=="C"' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	### subset variants with ALT == A, T, G, or C
	awk '$5=="A" || $5=="T" || $5=="G" || $5=="C"' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	# Split INFO column by ";"
	awk -F";" '$1=$1' OFS="\t" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	# Remove unwanted columns (new header: CHROM POS ID REF ALT QUAL FILTER AC AF EAS_AF AMR_AF AFR_AF EUR_AF SAS_AF  
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	# Remove unwanted strings
	awk '{ gsub("SAS_AF=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	awk '{ gsub("EUR_AF=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	awk '{ gsub("AFR_AF=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	awk '{ gsub("AMR_AF=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	awk '{ gsub("EAS_AF=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	awk '{ gsub("AF=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
	awk '{ gsub("AC=", ""); print }' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" && mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf

	# rename file 
	mv "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf "$OUT_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf

	echo "Chromosome $SGE_TASK_ID complete!"
	
