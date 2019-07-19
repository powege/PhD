#!/bin/bash

### Script that QCs gnomAD vcf:
# -- subsets SNVs with "PASS" filter status

### get controls_AC and controls_AF

# Set variables
RAW_ROOT=/well/lindgren/George/Data/gnomAD/vcf_raw
OUTPUT_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP

# Save # lines 
awk '/^#/' gnomad.genomes.r2.1.sites.chr21.vcf > gnomad.genomes.r2.1.sites.chr21.vcf.meta

# create dummy dataset
#tail -1000 gnomad.genomes.r2.1.sites.chr21.vcf > dummy.vcf

# grep from INFO colum into tmp files
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
}'  dummy.vcf > INFO_AF.tmp
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
}'  dummy.vcf > INFO_AF_female.tmp
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
}'  dummy.vcf > INFO_AF_female_control.tmp
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
}'  dummy.vcf > INFO_AF_female_non_neuro.tmp
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
}'  dummy.vcf > INFO_AF_male.tmp
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
}'  dummy.vcf > INFO_AF_male_control.tmp
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
}'  dummy.vcf > INFO_AF_male_non_neuro.tmp
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
}'  dummy.vcf > INFO_VQSLOD.tmp
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
}'  dummy.vcf > INFO_InbreedingCoeff.tmp
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
}'  dummy.vcf > INFO_n_alt_alleles.tmp

# remove text strings
awk '{ gsub("AF=", ""); print }' INFO_AF.tmp > tmp && mv tmp INFO_AF.tmp
awk '{ gsub("AF_male=", ""); print }' INFO_AF_male.tmp > tmp && mv tmp INFO_AF_male.tmp
awk '{ gsub("AF_female=", ""); print }' INFO_AF_female.tmp > tmp && mv tmp INFO_AF_female.tmp
awk '{ gsub("controls_AF_male=", ""); print }' INFO_AF_male_control.tmp > tmp && mv tmp INFO_AF_male_control.tmp
awk '{ gsub("controls_AF_female=", ""); print }' INFO_AF_female_control.tmp > tmp && mv tmp INFO_AF_female_control.tmp
awk '{ gsub("non_neuro_AF_male=", ""); print }' INFO_AF_male_non_neuro.tmp > tmp && mv tmp INFO_AF_male_non_neuro.tmp
awk '{ gsub("non_neuro_AF_female=", ""); print }' INFO_AF_female_non_neuro.tmp > tmp && mv tmp INFO_AF_female_non_neuro.tmp
awk '{ gsub("InbreedingCoeff=", ""); print }' INFO_InbreedingCoeff.tmp > tmp && mv tmp INFO_InbreedingCoeff.tmp
awk '{ gsub("VQSLOD=", ""); print }' INFO_VQSLOD.tmp > tmp && mv tmp INFO_VQSLOD.tmp
awk '{ gsub("n_alt_alleles=", ""); print }' INFO_n_alt_alleles.tmp > tmp && mv tmp INFO_n_alt_alleles.tmp

# cbind tmp INFO files
paste -d " " INFO_AF.tmp INFO_AF_male.tmp INFO_AF_female.tmp INFO_AF_male_non_neuro.tmp INFO_AF_female_non_neuro.tmp INFO_AF_male_control.tmp INFO_AF_female_control.tmp INFO_InbreedingCoeff.tmp INFO_VQSLOD.tmp INFO_n_alt_alleles.tmp > INFO.tmp

# Subset vcf columns: CHROM POS ID REF ALT QUAL FILTER
awk '{print $1, $2, $3, $4, $5, $6, $7}' dummy.vcf > dummy_formatted.vcf

# cbind INFO and vcf files
paste dummy_formatted.vcf INFO.tmp > tmp && mv tmp dummy_formatted.vcf

# remove tmp AF files 
rm INFO*
rm dummy.vcf

# Subset variants with "PASS" filter status (this removes header and meta lines)
awk '$7 == "PASS"' dummy_formatted.vcf > tmp && mv tmp dummy_formatted.vcf
#awk '$7 == "PASS" || $7 == "RF"' dummy_formatted.vcf > tmp && mv tmp dummy_formatted.vcf


