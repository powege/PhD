#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-22 -tc 16
#$ -N thousandGP_hom_prep_for_VEP
#$ -o /well/lindgren/George/Workflows/NC_constraint/Cluster_log/
#$ -e /well/lindgren/George/Workflows/NC_constraint/Cluster_log/

### Script that QCs 1000GP vcfs by chromosome.
### QC: PASS filter statuse; SNV; homozygous genotype >=1 
### Output: .vcf file (CHROM POS ID REF ALT QUAL FILTER AC AF EAS_AF AMR_AF AFR_AF EUR_AF SAS_AF)

# Set variables
RAW_ROOT=/well/lindgren/George/Data/1000GP/vcf_raw/
OUTPUT_ROOT=/well/lindgren/George/Data/1000GP/vcf_QCed_VEP/

#SGE_TASK_ID=_dummy
#head -10000 "$RAW_ROOT"ALL.chr17_GRCh38.genotypes.20170504.vcf > "$RAW_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38.genotypes.20170504.vcf

# Save # lines 
awk '/^#/' "$RAW_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38.genotypes.20170504.vcf > "$RAW_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38.genotypes.20170504.meta

# Remove meta data 
sed '/^#/ d' < "$RAW_ROOT"ALL.chr"$SGE_TASK_ID"_GRCh38.genotypes.20170504.vcf > "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
	#grep -o '^[^#]*' $FILE > "$FILE".QC_hom.temp
	
# Subset all rows with >= 1 individual with homozygous alternate allele "1|1"
awk '{ for (i=9;i<=NF;i++) if ($i == "1|1") { print; next } }'  "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID"
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
	
# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER INFO
awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID"
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf

# Subset variants with "PASS" filter status
awk '$7 == "PASS"' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp 
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf

# Subset variants with alternate base == 1
awk 'length($5)<=1' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp 
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf

# Subset variants with reference base == 1
awk 'length($4)<=1' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf

# Split INFO column by ";"
awk -F";" '$1=$1' OFS="\t" "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp 
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf

# Remove unwanted columns (new header: CHROM POS ID REF ALT QUAL FILTER AC AF EAS_AF AMR_AF AFR_AF EUR_AF SAS_AF  
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13, $14, $15, $16, $17}' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp 
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf

# Remove unwanted strings
awk '{ gsub("SAS_AF=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
awk '{ gsub("EUR_AF=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
awk '{ gsub("AFR_AF=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
awk '{ gsub("AMR_AF=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
awk '{ gsub("EAS_AF=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
awk '{ gsub("AF=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf
awk '{ gsub("AC=", ""); print }' "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf > tmp 
mv tmp "$RAW_ROOT"1000GP_VEP_input_hom_chr"$SGE_TASK_ID".vcf



