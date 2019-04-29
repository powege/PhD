#!/bin/bash

### Script that downloads and unzips MGP SNV vcf 

# set file directory
PED_ROOT=/well/lindgren/George/Data/MGP/vcf_raw/

# download MGP SNV vcf from ftp (ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz)
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
-P "$PED_ROOT"

# unzip file
gunzip "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf.gz

# save # lines 
awk '/^#/' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf > "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf.meta

# remove # lines
sed '/^#/ d' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf > "$PED_ROOT"tmp 
mv "$PED_ROOT"tmp "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf

# split by chromosome
awk '{print $0>$1 "_mgp.v5.merged.snps_all.dbSNP142.vcf"}' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
mv "$PED_ROOT""$CHR"_mgp.v5.merged.snps_all.dbSNP142.vcf "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142_chr"$CHR".vcf
done

# rm X, MT, and Y chr
rm X_mgp.v5.merged.snps_all.dbSNP142.vcf 
rm Y_mgp.v5.merged.snps_all.dbSNP142.vcf 
rm MT_mgp.v5.merged.snps_all.dbSNP142.vcf 




