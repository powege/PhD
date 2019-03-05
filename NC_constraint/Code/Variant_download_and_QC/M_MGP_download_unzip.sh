#!/bin/bash

### Script that downloads and unzips MGP SNV vcf 

# Set working directory for files to download into
cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/Mouse/MGP

# Download MGP SNV vcf from ftp (ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz)
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz

# unzip file
gunzip mgp.v5.merged.snps_all.dbSNP142.vcf.gz


