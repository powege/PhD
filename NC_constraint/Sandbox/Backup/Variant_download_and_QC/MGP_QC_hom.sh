#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -t 1-19 -tc 16
#$ -N MGP_vcf_QC
#$ -o /well/lindgren/George/Workflows/NC_constraint/Log/
#$ -e /well/lindgren/George/Workflows/NC_constraint/Log/

### Script that QCs MGP vcf   -- subsets SNV sites with one or more high confidence homozygous genotype call 
#                                 (denoted with an FI tag of 1)

# Set directory
PED_ROOT=/well/lindgren/George/Data/MGP/vcf_raw/
OUT_ROOT=/well/lindgren/George/Data/MGP/vcf_QCed_VEP/

# SGE_TASK_ID=19
# head -10000 "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142_chr1.vcf > "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142_chr"$SGE_TASK_ID".vcf

# Subset variants with alternate base == 1
awk 'length($5)<=1' "$PED_ROOT"mgp.v5.merged.snps_all.dbSNP142_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" 
mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf

# Subset variants with reference base == 1
awk 'length($4)<=1' "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" 
mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf

# subset SNV sites with one or more high confidence homozygous genotype call 
grep -P -a "^.+?\s1\/1.+" "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" 
mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf

# Subset vcf columns: CHROM POS ID REF ALT FILTER (ie remove INFO)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"tmp_"$SGE_TASK_ID" 
mv "$PED_ROOT"tmp_"$SGE_TASK_ID" "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf
mv "$PED_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf "$OUT_ROOT"MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf

#####

# for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
# do
# wc -l MGP_v5_snps_QCed_hom_chr"$CHR".vcf
# done

### raw colnames:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	129P2_OlaHsd	129S1_SvImJ	129S5SvEvBrd	AKR_J	A_J	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J	CAST_EiJ	CBA_J	DBA_1J	DBA_2J	FVB_NJ	I_LnJ	KK_HiJ	LEWES_EiJ	LP_J	MOLF_EiJ	NOD_ShiLtJ	NZB_B1NJ	NZO_HlLtJ	NZW_LacJ	PWK_PhJ	RF_J	SEA_GnJ	SPRET_EiJ	ST_bJ	WSB_EiJ	ZALENDE_EiJ




