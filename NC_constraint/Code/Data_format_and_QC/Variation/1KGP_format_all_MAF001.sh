#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -t 1-22 -tc 16
#$ -N thouGP_MAF001
#$ -o /well/lindgren/George/Workflows/NC_constraint/Log/
#$ -e /well/lindgren/George/Workflows/NC_constraint/Log/

### Script that QCs 1000GP vcfs by chromosome.
### QC: PASS filter status; SNV 
### Output cols: 1.CHROM 2.POS 3.ID 4.REF 5.ALT 6.QUAL 7.FILTER 8.AC 9.AF 

### Set working directory
PED_ROOT=/well/lindgren/George/Data/1KGP/Variants/vcf_QCed_VEP/

# SGE_TASK_ID=dummy
# head -10000 "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr17.vcf > "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf

	echo "Processing Chromosome $SGE_TASK_ID"

# Subset SNVs with MAF >= 0.001
awk '$9 >= 0.001' "$PED_ROOT"1000GP_phase3_snvs_QCed_all_chr"$SGE_TASK_ID".vcf > "$PED_ROOT"1000GP_phase3_snvs_QCed_all_MAF001_chr"$SGE_TASK_ID".vcf
	
	echo "Chromosome $SGE_TASK_ID complete!"
	
	
	
	
	
	
	