#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -N wild_allSPECIES_vcf_QC
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that QCs wild mouse vcf  

# Set directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/Wild_mice/Variants/Raw/
OUT_ROOT=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/Wild_mice/Variants/Formatted/
    
# Set outfilenames
IN_FILE=AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf
TMP_FILE=tmp_file

# test file
# head -10000 "$PED_ROOT"AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf > "$PED_ROOT"test_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf
    
# unzip 
gunzip "$PED_ROOT""$IN_FILE".gz
    
# Remove # lines
awk '!/^#/' "$PED_ROOT""$IN_FILE" > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT""$TMP_FILE"

# Subset vcf columns for output: CHROM POS ID REF ALT FILTER 
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$PED_ROOT""$TMP_FILE" > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT""$TMP_FILE"
    
# Subset variants with alternate base == 1
awk 'length($5)<=1' "$PED_ROOT""$TMP_FILE" > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT""$TMP_FILE"
    
# Subset variants with reference base == 1
awk 'length($4)<=1' "$PED_ROOT""$TMP_FILE" > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT""$TMP_FILE"

# Split by chromosome
awk '{print $0>$1 "_wild_mice_all_snps_QCed.vcf"}' "$PED_ROOT""$TMP_FILE"
for CHR in {1..19} X
do
mv "$PED_ROOT"chr"$CHR"_wild_mice_all_snps_QCed.vcf "$OUT_ROOT"Wild_mice_allSPECIES_snps_PASS_chr"$CHR".vcf
done

# Zip 
for CHR in {1..19} X
do
gzip "$OUT_ROOT"Wild_mice_allSPECIES_snps_PASS_chr"$CHR".vcf
done

# Zip 
gzip "$PED_ROOT""$IN_FILE"

# Remove tmp files 
rm "$PED_ROOT""$TMP_FILE"
rm "$PED_ROOT"chr*_wild_mice_all_snps_QCed.vcf




    