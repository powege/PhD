#!/bin/bash

### Script that QCs 1000GP vcfs by chromosome.
### QC: PASS filter statuse; SNV 
### Output cols: 1.CHROM 2.POS 3.ID 4.REF 5.ALT 6.QUAL 7.FILTER 8.AC 9.AF 10.EAS_AF 11.AMR_AF 12.AFR_AF 13.EUR_AF 14.SAS_AF

### Set working directory
PED_ROOT=/NGS/users/George/Workflows/PC_constraint/Paper/Data/1000GP/
#cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/Human/1000GP/


### FOR LOOP

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do 
	echo "Processing Chromosome $CHR"

	# Save # lines 
	awk '/^#/' "$PED_ROOT"ALL.chr"$CHR"_GRCh38_sites.20170504.vcf  > "$PED_ROOT"ALL.chr"$CHR"_GRCh38_sites.20170504.vcf.meta
	# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER INFO
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' "$PED_ROOT"ALL.chr"$CHR"_GRCh38_sites.20170504.vcf > "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	# Subset variants with "PASS" filter status (this removes header and meta lines)
	awk '$7 == "PASS"' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	# Subset variants with alternate base == 1
	awk 'length($5)<=1' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	# Subset variants with reference base == 1
	awk 'length($4)<=1' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	### subset variants with REF == A, T, G, or C
	awk '$4=="A" || $4=="T" || $4=="G" || $4=="C"' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	### subset variants with ALT == A, T, G, or C
	awk '$5=="A" || $5=="T" || $5=="G" || $5=="C"' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	# Split INFO column by ";"
	awk -F";" '$1=$1' OFS="\t" "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	# Remove unwanted columns (new header: CHROM POS ID REF ALT QUAL FILTER AC AF EAS_AF AMR_AF AFR_AF EUR_AF SAS_AF  
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $14, $15, $16, $17}' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	# Remove unwanted strings
	awk '{ gsub("SAS_AF=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	awk '{ gsub("EUR_AF=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	awk '{ gsub("AFR_AF=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	awk '{ gsub("AMR_AF=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	awk '{ gsub("EAS_AF=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	awk '{ gsub("AF=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf
	awk '{ gsub("AC=", ""); print }' "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf

	# rename file 
	"$PED_ROOT"H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf >> "$PED_ROOT"H_1000GP_QCed_VEP_input_all.vcf
	
	echo "Chromosome $CHR complete!"
done
