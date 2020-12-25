#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-19 -tc 16
#$ -N MGP_format
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/


echo "########################################################"
echo "Submit sample QC job"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

### LOAD MODULE
module load R/3.4.3

### SET VARS 

R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Thesis_code/Download_format_QC/MGP/Variants/2.MGP_vcf_format_QC.R

IN_VCF=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/MGP/Variants/Raw/mgp.v5.merged.snps_all.dbSNP142_chr"$SGE_TASK_ID".vcf.gz
OUT_VCF_ALL=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/MGP/Variants/Formatted/MGP_v5_all_snps_PASS_chr"$SGE_TASK_ID".vcf.gz
OUT_VCF_CLASSICAL=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/MGP/Variants/Formatted/MGP_v5_classical_snps_PASS_chr"$SGE_TASK_ID".vcf.gz
OUT_VCF_WILD=/well/lindgren/George/Data/Thesis_workflow/Thesis_data/MGP/Variants/Formatted/MGP_v5_wild_snps_PASS_chr"$SGE_TASK_ID".vcf.gz

### RUN

Rscript --vanilla $R_SCRIPT \
$IN_VCF \
$OUT_VCF_ALL \
$OUT_VCF_CLASSICAL \
$OUT_VCF_WILD

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0