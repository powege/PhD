#!/bin/bash

### Script that downloads MGP SNV vcf from FTP
# splits into gzipped chromosome-specifc files

# set file directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/MGP/Variants/Raw/

# download MGP SNV vcf from ftp
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
-P "$PED_ROOT"

# unzip file
gunzip -c "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf.gz > "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf

# save # lines
awk '/^#/' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf > "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf.meta

# remove # lines
sed '/^#/ d' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf

# split by chromosome
cd $PED_ROOT
awk '{print $0>$1 "_mgp.v5.merged.snps_all.dbSNP142.vcf"}' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf
for CHR in {1..19} X
do
mv "$PED_ROOT""$CHR"_mgp.v5.merged.snps_all.dbSNP142.vcf "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142_chr"$CHR".vcf
done

# gzip files
for CHR in {1..19} X
do
gzip "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142_chr"$CHR".vcf
echo "$CHR"
done

# rm MT and Y chr and unzipped all
rm "$PED_ROOT"*_mgp.v5.merged.snps_all.dbSNP142.vcf 
rm "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142.vcf



