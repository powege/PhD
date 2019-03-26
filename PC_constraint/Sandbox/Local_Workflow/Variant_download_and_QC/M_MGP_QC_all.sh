#!/bin/bash

### Script that QCs MGP vcf -- subsets SNVs with "PASS" filter status

# Set working directory
cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/Mouse/MGP/

# Save # lines 
awk '/^#/' mgp.v5.merged.snps_all.dbSNP142.vcf > mgp.v5.merged.snps_all.dbSNP142.vcf.meta

# Subset vcf columns: CHROM POS ID REF ALT FILTER (ie remove INFO)
awk '{print $1, $2, $3, $4, $5, $6, $7}' mgp.v5.merged.snps_all.dbSNP142.vcf > MGP_temp

# Subset variants with "PASS" filter status (this removes header and meta lines)
awk '$7 == "PASS"' MGP_temp > tmp && mv tmp MGP_temp

# Subset vcf columns: CHROM POS ID REF ALT
awk '{print $1, $2, $3, $4, $5, $6, $7}' MGP_temp > tmp && mv tmp MGP_temp

# Subset variants with alternate base == 1
awk 'length($5)<=1' MGP_temp > tmp && mv tmp MGP_temp

# Subset variants with reference base == 1
awk 'length($4)<=1' MGP_temp > tmp && mv tmp MGP_temp

# Subset autosomes and X chromosome
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" || $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  || $1 == "17" || $1 == "18" || $1 == "19"  || $1 == "X")' MGP_temp > tmp && mv tmp MGP_temp

# Rename output
mv MGP_temp M_MGP_QCed_VEP_input_all.vcf

# Split vcf by chromosome
awk '{print $0>$1 "_M_chr_QCed.vcf"}' M_MGP_QCed_all.vcf
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X
do
mv "$CHR"_M_chr_QCed.vcf M_MGP_QCed_VEP_input_all_chr"$CHR".vcf
done




