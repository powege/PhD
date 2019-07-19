#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-22 -tc 16
#$ -N gnomAD_VEP_run
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that runs Ensembl VEP v94 on formatted .vcf files

### load the VEP module
module load ensembl-tools/94

### set rootnames of your files
IN_ROOT=/well/lindgren/George/Data/gnomAD/vcf_raw/
OUT_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/

### set file names
IN_FILE=gnomAD_v2.1.1_GRC38_snps_QCed_controls_chr
OUT_FILE=gnomAD_v2.1.1_GRC38_snps_QCed_controls_VEPout_chr
  
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
$VEP --input_file "$IN_ROOT""$IN_FILE""$SGE_TASK_ID".vcf  \
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
--output_file "$OUT_ROOT""$OUT_FILE""$SGE_TASK_ID".vcf \
--no_stats \
--offline

