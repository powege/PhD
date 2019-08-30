#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-19 -tc 16
#$ -N MGP_VEP_run
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that runs Ensembl VEP v94 on formatted .vcf files

### load the VEP module
module load ensembl-tools/94

### set rootnames of your files
IN_ROOT=/well/lindgren/George/Data/MGP/Variants/vcf_QCed_VEP/
OUT_ROOT=/well/lindgren/George/Data/MGP/Variants/vcf_QCed_VEP/

### set file names
IN_FILE=MGP_v5_allSTRAIN_snps_QCed_hom_chr
OUT_FILE=MGP_v5_allSTRAIN_snps_QCed_VEPout_hom_chr
# IN_FILE=MGP_v5_allMUSMUS_snps_QCed_hom_chr
# OUT_FILE=MGP_v5_allMUSMUS_snps_QCed_VEPout_hom_chr
  
  
  ### set your VEP variable
  VEP=/well/lindgren/George/ensembl-vep/vep

### make sure VEP version 94 is git cloned
#cd /well/lindgren/George/
#git clone https://github.com/Ensembl/ensembl-vep
#cd ensembl-vep
#git checkout release/94
### Create .vep dir and run INSTALL.pl to download cache files
#mkdir /well/lindgren/George/.vep
# perl INSTALL.pl \
# --AUTO c \
# --ASSEMBLY "GRCh38" \
# --CACHEDIR /well/lindgren/George/.vep/ \
# --SPECIES mus_musculus \
# --VERSION 94
### or wget and unzip cache files from ftp into .vep:
# https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
# ftp://ftp.ensembl.org/pub/release-94/variation/indexed_vep_cache/

# check files are in place
#ls $HOME/.vep/
#ls /well/lindgren/George/.vep/

### RUN VEP
# -- input file:   "$IN_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf
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
--species "mus_musculus" \
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

