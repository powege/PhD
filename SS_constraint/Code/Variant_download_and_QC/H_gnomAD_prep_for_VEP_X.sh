#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -N gnomAD_prep_for_VEP_X
#$ -o /well/lindgren/George/Workflows/SS_constraint/Log/
#$ -e /well/lindgren/George/Workflows/SS_constraint/Log/


### Script that preps gnomAD vcf for VEP:
# -- extracts from INFO column:
# - AF; 
# - AF_male; 
# - AF_female; 
# - non_neuro_AF_male; 
# - non_neuro_AF_female; 
# - controls_AF_male; 
# - controls_AF_female; 
# - InbreedingCoeff; 
# - VQSLOD; 
# - n_alt_alleles
# -- subsets variants with "PASS" filter status

# Set variables
RAW_ROOT=/well/lindgren/George/Data/gnomAD/vcf_raw/
OUTPUT_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/
SGE_TASK_ID=X

# create dummy dataset
#head -10000 gnomad.genomes.r2.1.sites.chr21.vcf > dummy.vcf

# Save # lines 
#awk '/^#/' "$RAW_ROOT"gnomad.genomes.r2.1.sites.chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"gnomad.genomes.r2.1.sites.chr"$SGE_TASK_ID".vcf.meta

# Remove # lines
awk '!/^#/' "$RAW_ROOT"gnomad.genomes.r2.1.sites.chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf

# awk from INFO colum into tmp files
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^AF=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^AF_female=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_female_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^controls_AF_female=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_female_control_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^non_neuro_AF_female=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_female_non_neuro_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^AF_male=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_male_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^controls_AF_male=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_male_control_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^non_neuro_AF_male=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_male_non_neuro_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^VQSLOD=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^InbreedingCoeff=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^n_alt_alleles=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp

# remove text strings
awk '{ gsub("AF=", ""); print }' "$RAW_ROOT"INFO_AF_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \ 
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_chr"$SGE_TASK_ID".tmp
awk '{ gsub("AF_male=", ""); print }' "$RAW_ROOT"INFO_AF_male_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_male_chr"$SGE_TASK_ID".tmp
awk '{ gsub("AF_female=", ""); print }' "$RAW_ROOT"INFO_AF_female_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_female_chr"$SGE_TASK_ID".tmp
awk '{ gsub("controls_AF_male=", ""); print }' "$RAW_ROOT"INFO_AF_male_control_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_male_control_chr"$SGE_TASK_ID".tmp
awk '{ gsub("controls_AF_female=", ""); print }' "$RAW_ROOT"INFO_AF_female_control_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_female_control_chr"$SGE_TASK_ID".tmp
awk '{ gsub("non_neuro_AF_male=", ""); print }' "$RAW_ROOT"INFO_AF_male_non_neuro_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_male_non_neuro_chr"$SGE_TASK_ID".tmp
awk '{ gsub("non_neuro_AF_female=", ""); print }' "$RAW_ROOT"INFO_AF_female_non_neuro_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_female_non_neuro_chr"$SGE_TASK_ID".tmp
awk '{ gsub("InbreedingCoeff=", ""); print }' "$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp
awk '{ gsub("VQSLOD=", ""); print }' "$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp
awk '{ gsub("n_alt_alleles=", ""); print }' "$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp

# cbind tmp INFO files
paste -d " " "$RAW_ROOT"INFO_AF_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_male_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_female_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_male_non_neuro_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_female_non_neuro_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_male_control_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_female_control_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"INFO_chr"$SGE_TASK_ID".tmp

# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf

# cbind INFO and vcf files
paste -d " " "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf "$RAW_ROOT"INFO_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf

# remove tmp files 
rm "$RAW_ROOT"*_chr"$SGE_TASK_ID".tmp

# Subset variants with "PASS" filter status (this removes header and meta lines)
awk '$7 == "PASS"' "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf
#awk '$7 == "PASS" || $7 == "RF"' dummy_formatted.vcf > tmp && mv tmp dummy_formatted.vcf


