#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -t 1-19 -tc 16
#$ -N MGP_all_strain_n_DP
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

R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/MGP/Coverage/4.MGP_all_strain_n_DP_less10X.R

CLASSICAL_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Formatted/MGP_classical_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed
WILD_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Formatted/MGP_wild_derived_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R
CHR="$SGE_TASK_ID"
OUT_BED_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Formatted/MGP_all_strain_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed

### RUN

Rscript --vanilla $R_SCRIPT \
$CLASSICAL_FILE \
$WILD_FILE \
$FUNCTIONS_FILE \
$CHR \
$OUT_BED_FILE 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0