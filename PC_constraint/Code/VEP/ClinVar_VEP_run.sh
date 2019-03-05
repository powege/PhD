#!/bin/bash

### Script that runs Ensembl VEP on formatted ClinVar vcf

### cd into vcf directory
cd /Users/g.powell/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP

### $PED_ROOT is the rootname of your files
PED_ROOT=/Users/g.powell/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP

### Run the perl INSTALL.pl script contained in the docker image
#docker run -t -i -v $PED_ROOT:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl

### Run the docker image and load bash
docker run -t -i -v $PED_ROOT:/opt/vep/.vep ensemblorg/ensembl-vep /bin/bash

# check files are in place
ls /opt/vep/.vep/

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
./vep -i /opt/vep/.vep/ClinVar_VEP_input.vcf --cache --ccds --symbol --biotype --canonical --vcf --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CCDS" --pick -o /opt/vep/.vep/ClinVar_VEP_output.vcf --stats_file /opt/vep/.vep/ClinVar_VEP_output.html --offline
