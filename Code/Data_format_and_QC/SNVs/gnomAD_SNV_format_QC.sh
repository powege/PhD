#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-22 -tc 16
#$ -N gnomAD_controls_QC
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/


### Script that preps gnomAD vcf for VEP:
# -- extracts from INFO column:
# - controls_AC; 
# - controls_AF; 
# - InbreedingCoeff; 
# - VQSLOD; 
# - n_alt_alleles
# -- subsets SNVs with "PASS" filter status and MAF > 0
### Output columns:
#1 CHR
#2 POS
#3 ID 
#4 REF
#5 ALT
#6 QUAL
#7 FILTER
#8 controls_AC
#9 controls_AF
#10 InbreedingCoeff
#11 VQSLOD
#12 n_alt_alleles

# Set variables
RAW_ROOT=/well/lindgren/George/Data/gnomAD/vcf_raw/
OUTPUT_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/
OUT_FILE_NAME=gnomAD_v2.1.1_GRC38_snps_QCed_controls_chr

# create dummy dataset
# SGE_TASK_ID=dummy
# head -10000 gnomad.genomes.r2.1.sites.chr21.vcf > gnomad.genomes.r2.1.sites.chr"$SGE_TASK_ID".vcf

# Save # lines 
#awk '/^#/' "$RAW_ROOT"gnomad.genomes.r2.1.sites.chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"gnomad.genomes.r2.1.sites.chr"$SGE_TASK_ID".vcf.meta

# Remove # lines
awk '!/^#/' "$RAW_ROOT"gnomad.genomes.r2.1.1.sites."$SGE_TASK_ID".liftover_grch38.vcf > "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf

# awk from INFO colum into tmp files
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^controls_AC=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AC_control_chr"$SGE_TASK_ID".tmp
awk -F';' '
{
  val=""
  for(i=1;i<=NF;i++){
     if($i ~ /^controls_AF=[0-9]+/){
         val=(val?val OFS $i:$i)
     }
  }
  if(val){
     print val
  }
  else{
     print "NA"
  }
}'  "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_AF_control_chr"$SGE_TASK_ID".tmp
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
}'  "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp
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
}'  "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp
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
}'  "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp

# remove text strings
awk '{ gsub("controls_AC=", ""); print }' "$RAW_ROOT"INFO_AC_control_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AC_control_chr"$SGE_TASK_ID".tmp
awk '{ gsub("controls_AF=", ""); print }' "$RAW_ROOT"INFO_AF_control_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_AF_control_chr"$SGE_TASK_ID".tmp
awk '{ gsub("InbreedingCoeff=", ""); print }' "$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp
awk '{ gsub("VQSLOD=", ""); print }' "$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp
awk '{ gsub("n_alt_alleles=", ""); print }' "$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp

# cbind tmp INFO files
paste -d " " "$RAW_ROOT"INFO_AC_control_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_AF_control_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_InbreedingCoeff_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_VQSLOD_chr"$SGE_TASK_ID".tmp \
"$RAW_ROOT"INFO_n_alt_alleles_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"INFO_chr"$SGE_TASK_ID".tmp

# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf

# cbind INFO and vcf files
paste -d " " "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf "$RAW_ROOT"INFO_chr"$SGE_TASK_ID".tmp > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf

# remove tmp files 
rm "$RAW_ROOT"*_chr"$SGE_TASK_ID".tmp

# Subset variants with "PASS" filter status (this removes header and meta lines)
awk '$7 == "PASS"' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf
#awk '$7 == "PASS" || $7 == "RF"' dummy_formatted.vcf > tmp && mv tmp dummy_formatted.vcf

# Subset variants with alternate base == 1
awk 'length($5)<=1' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf
# Subset variants with reference base == 1
awk 'length($4)<=1' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf
### subset variants with REF == A, T, G, or C
awk '$4=="A" || $4=="T" || $4=="G" || $4=="C"' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf
### subset variants with ALT == A, T, G, or C
awk '$5=="A" || $5=="T" || $5=="G" || $5=="C"' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf
### subset variants with MAF > 0 
awk '$9>0' "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf > "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" && \
mv "$RAW_ROOT"tmp_chr"$SGE_TASK_ID" "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf


# rename output
mv "$RAW_ROOT"gnomAD_tmp_chr"$SGE_TASK_ID".vcf "$RAW_ROOT""$OUT_FILE_NAME""$SGE_TASK_ID".vcf

