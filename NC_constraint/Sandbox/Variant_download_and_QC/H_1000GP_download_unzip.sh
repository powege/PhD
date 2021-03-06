#!/bin/bash

### Script that downloads and unzips 1000GP vcf 

# set file directory
PED_ROOT=/well/lindgren/George/Data/1000GP/vcf_raw/

# Download 1000GP SNP vcf from ftp (GRCh38)
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr"$CHR"_GRCh38_sites.20170504.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr"$CHR"_GRCh38.genotypes.20170504.vcf.gz \
-P "$PED_ROOT"
done

# unzip files
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
#gunzip ALL.chr"$CHR"_GRCh38_sites.20170504.vcf.gz
gunzip "$PED_ROOT"ALL.chr"$CHR"_GRCh38.genotypes.20170504.vcf.gz
done


### OLD FILES IN SAND
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#gunzip ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


