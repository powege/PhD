#!/bin/bash

### Script that QCs 1000GP vcfs by chromosome.
### QC: PASS filter statuse; SNV; homozygous genotype >=1 
### Output: .vcf file (CHROM POS ID REF ALT QUAL FILTER); .af file (AC AF EAS_AF AMR_AF AFR_AF EUR_AF SAS_AF)

### Set working directory
cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/Human/1000GP/

### FOR LOOP

for FILE in *.genotypes.vcf
do
	echo "Processing $FILE"
	
	# Save # lines 
	awk '/^#/' $FILE > "$FILE".meta
	# Remove meta data 
	sed '/^#/ d' < $FILE > "$FILE".QC_hom.temp
	#grep -o '^[^#]*' $FILE > "$FILE".QC_hom.temp
	# Subset all rows with >= 1 individual with homozygous alternate allele "1|1"
	awk '{ for (i=9;i<=NF;i++) if ($i == "1|1") { print; next } }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER INFO
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	# Subset variants with "PASS" filter status
	awk '$7 == "PASS"' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	# Subset variants with alternate base == 1
	awk 'length($5)<=1' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	# Subset variants with reference base == 1
	awk 'length($4)<=1' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp

	# Split INFO column by ";"
	awk -F";" '$1=$1' OFS="\t" "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	# Remove unwanted columns (new header: CHROM POS ID REF ALT QUAL FILTER AC AF EAS_AF AMR_AF AFR_AF EUR_AF SAS_AF  
	awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $14, $15, $16, $17}' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	# Remove unwanted strings
	awk '{ gsub("SAS_AF=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	awk '{ gsub("EUR_AF=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	awk '{ gsub("AFR_AF=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	awk '{ gsub("AMR_AF=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	awk '{ gsub("EAS_AF=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	awk '{ gsub("AF=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp
	awk '{ gsub("AC=", ""); print }' "$FILE".QC_hom.temp > tmp && mv tmp "$FILE".QC_hom.temp

	# Split into vcf and af files
	#awk '{print $1, $2, $3, $4, $5, $6, $7}' "$FILE".QC_hom.temp >> H_1000GP_QCed_hom.vcf
	#awk '{print $8, $9, $10, $11, $12, $13, $14}' "$FILE".QC_hom.temp >> H_1000GP_QCed_hom.af
	
	# rename file 
    "$FILE".QC_hom.temp > H_1000GP_QCed_VEP_input_hom.vcf
	
	# remove temp file
	rm "$FILE".QC_hom.temp 
	
	echo "$FILE complete!"
done

