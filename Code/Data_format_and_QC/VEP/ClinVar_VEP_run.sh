#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -N ClinVar_VEP_run
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that runs Ensembl VEP v94 on formatted .vcf files

### load the VEP module
module load ensembl-tools/94

### set rootnames of your files
IN_ROOT=/well/lindgren/George/Data/ClinVar/formatted/
OUT_ROOT=/well/lindgren/George/Data/ClinVar/formatted/

### set file names
IN_FILE=ClinVar_pathogenic_formatted.vcf
OUT_FILE=ClinVar_pathogenic_snps_QCed_VEPout.vcf
  
  ### set your VEP variable
  VEP=/well/lindgren/George/ensembl-vep/vep

### make sure VEP version 94 is git cloned
#cd /well/lindgren/George/
#git clone https://github.com/Ensembl/ensembl-vep
#cd ensembl-vep
#git checkout release/94
#mkdir /well/lindgren/George/.vep
#perl INSTALL.pl \
#--AUTO c \
#--ASSEMBLY "GRCh38" \
#--CACHEDIR /well/lindgren/George/.vep/ \
#--SPECIES "homo_sapien" \
#--VERSION 94

# check files are in place
#ls $HOME/.vep/
#ls /well/lindgren/George/.vep/

### RUN VEP
# -i -- input file:   1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf
# Flags: (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#output)
# --cache --
# --dir_cache --
# --ccds --
# --symbol -- 
# --biotype --
# --canonical --
# --vcf -- output as vcf file
# --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" -- specify output fields
# --pick -- one annotation per variant (https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick)
$VEP --input_file "$IN_ROOT""$IN_FILE" \
--dir_cache /well/lindgren/George/.vep \
--cache \
--species "homo_sapiens" \
--ccds \
--symbol \
--biotype \
--canonical \
--vcf \
--fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" \
--pick \
--output_file "$OUT_ROOT""$OUT_FILE" \
--no_stats \
--offline

