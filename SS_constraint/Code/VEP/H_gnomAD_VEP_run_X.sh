#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -N gnomAD_VEP_run_X
#$ -o /well/lindgren/George/Workflows/SS_constraint/Log/
#$ -e /well/lindgren/George/Workflows/SS_constraint/Log/

### Script that runs Ensembl VEP on formatted gnomAD_formatted_for_VEP_chr"$CHR".vcf files

### load the VEP module
module load ensembl-tools/94

### $PED_ROOT is the rootname of your files
IN_ROOT=/well/lindgren/George/Data/gnomAD/vcf_raw/
OUT_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/
SGE_TASK_ID=X

### set your VEP variable
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
# -i -- input file: H_1000GP_QCed_dummy.vcf
# Flags: (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#output)
# --cache --
# --dir_cache --
# --ccds --
# --symbol -- 
# --biotype --
# --canonical --
# --vcf -- output as vcf file
# --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" -- specify output fields
# --pick -- one annotation per variant (https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick)
$VEP -i "$IN_ROOT"gnomAD_formatted_for_VEP_chr"$SGE_TASK_ID".vcf  \
 --dir_cache /well/lindgren/George/.vep \
 --cache \
 --ccds \
 --symbol \
 --biotype \
 --canonical \
 --vcf \
 --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" \
 --pick \
 -o "$OUT_ROOT"gnomAD_VEP_output_chr"$SGE_TASK_ID".vcf \
 --no_stats \
 --offline




#######



 ##--format "vcf" \


##for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
##do 
##./vep -i /opt/vep/.vep/H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf --cache --ccds --symbol --biotype --canonical --vcf --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" --pick -o /opt/vep/.vep/H_1000GP_QCed_VEP_output_all_chr"$CHR".vcf --stats_file /opt/vep/.vep/H_1000GP_QCed_VEP_output_all_chr"$CHR".html --offline
##done


## --compress_output
## --vcf -- output as vcf file


#### Run the perl INSTALL.pl script contained in the docker image
##docker run -t -i -v $PED_ROOT:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl

#### Run the docker image and load bash
##docker run -t -i -v $PED_ROOT:/opt/vep/.vep ensemblorg/ensembl-vep /bin/bash

#tail -n +1000000 gnomAD_formatted_for_VEP_chr17.vcf | head -n 10000 - > dummy_formatted_for_VEP_chr17.vcf 
#awk '{print $1, $2, $3, $4, $5}' dummy_formatted_for_VEP_chr17.vcf  > tmp && mv tmp dummy_formatted_for_VEP_chr17.vcf

##awk 'gsub(" ",",")' dummy_formatted_for_VEP_chr17.vcf > tmp
#awk '$5=="A" || $5=="T" || $5=="G" || $5=="C"' dummy_formatted_for_VEP_chr17.vcf > tmp
#awk '$4=="A" || $4=="T" || $4=="G" || $4=="C"' tmp > tmp2 && mv tmp2 tmp

