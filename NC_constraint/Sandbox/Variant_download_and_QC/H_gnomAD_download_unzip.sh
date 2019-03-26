#!/bin/bash

### Script that downloads and unzips gnomAD vcf 

# Set working directory for files to download into
#cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/gnomAD
cd /well/lindgren/George/Data/gnomAD/vcf_raw

# Download gnomAD genomes vcf from ftp
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr"$CHR".vcf.bgz
done

# unzip files
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
mv gnomad.genomes.r2.1.sites.chr"$CHR".vcf.bgz gnomad.genomes.r2.1.sites.chr"$CHR".vcf.gz
gunzip gnomad.genomes.r2.1.sites.chr"$CHR".vcf.gz
done

### OLD FILES IN SAND
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#gunzip ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


