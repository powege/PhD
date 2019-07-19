#!/bin/bash

### Script that downloads and unzips gnomAD vcf 

# Set working directory for files to download into
PED_ROOT=/well/lindgren/George/Data/gnomAD/vcf_raw/

# Download 1000GP SNP vcf from ftp (GRCh38)
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites."$CHR".liftover_grch38.vcf.bgz \
-P "$PED_ROOT"
done

# unzip files
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
mv "$PED_ROOT"gnomad.genomes.r2.1.1.sites."$CHR".liftover_grch38.vcf.bgz "$PED_ROOT"gnomad.genomes.r2.1.1.sites."$CHR".liftover_grch38.vcf.gz
gunzip "$PED_ROOT"gnomad.genomes.r2.1.1.sites."$CHR".liftover_grch38.vcf.gz
done

