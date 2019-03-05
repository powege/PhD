#!/bin/bash

### Script that runs Ensembl VEP on formatted gnomAD_formatted_for_VEP_chr"$CHR".vcf files

cd /apps/well/ensembl-tools/94/ensembl-vep-release-94/
perl INSTALL.pl -c /well/lindgren/George/

### cd into vcf directory
cd /well/lindgren/George/Data/gnomAD/vcf_QCed_VEP/

### $PED_ROOT is the rootname of your files
PED_ROOT=/well/lindgren/George/Data/gnomAD/vcf_QCed_VEP

# load the VEP module
module load ensembl-tools/94

# set your VEP variable
VEP=/apps/well/ensembl-tools/94/ensembl-vep-release-94/vep

### Run the perl INSTALL.pl script contained in the docker image
#docker run -t -i -v $PED_ROOT:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl

### Run the docker image and load bash
#docker run -t -i -v $PED_ROOT:/opt/vep/.vep ensemblorg/ensembl-vep /bin/bash

# check files are in place
#ls /opt/vep/.vep/

### RUN VEP
# -i -- input file: H_1000GP_QCed_dummy.vcf
# Flags: (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#output)
# --cache --
# --ccds --
# --symbol -- 
# --biotype --
# --canonical --
# --vcf -- output as vcf file
# --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" -- specify output fields
# --pick -- one annotation per variant (https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick)
$VEP -i dummy_formatted_for_VEP_chr17.vcf --cache --ccds --symbol --biotype --canonical --vcf --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" --pick -o dummy_VEP_output_chr17.vcf --stats_file dummy_VEP_output_chr17.html --offline


#for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
#do 
#./vep -i /opt/vep/.vep/H_1000GP_QCed_VEP_input_all_chr"$CHR".vcf --cache --ccds --symbol --biotype --canonical --vcf --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" --pick -o /opt/vep/.vep/H_1000GP_QCed_VEP_output_all_chr"$CHR".vcf --stats_file /opt/vep/.vep/H_1000GP_QCed_VEP_output_all_chr"$CHR".html --offline
#done


# --compress_output
# --vcf -- output as vcf file

