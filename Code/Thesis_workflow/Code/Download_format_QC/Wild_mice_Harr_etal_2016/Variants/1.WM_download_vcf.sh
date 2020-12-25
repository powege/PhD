#!/bin/bash

# Script that downloads vcf files from ftp 

# Set working directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/

# download all species with PASS vcf
wget http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz \
-P "$PED_ROOT" 

# download by population
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmc_CAST.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmd_FRA.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmd_GER.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmd_HEL.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmd_IRA.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmm_AFG.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmm_CZE.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Mmm_KAZ.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/pop_vcf/Ms_SPRE.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf.gz




