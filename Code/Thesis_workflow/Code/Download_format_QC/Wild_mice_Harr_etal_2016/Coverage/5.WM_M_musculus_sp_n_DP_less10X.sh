#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 4
#$ -t 1-19 -tc 16
#$ -N WM_musculus_sp_n_DP
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

R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Wild_mice_Harr_etal_2016/Coverage/5.WM_M_musculus_sp_n_DP_less10X.R

CASTANEUS_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_castaneus_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed.gz
DOMESTICUS_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_domesticus_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed.gz
MUSCULUS_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_musculus_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed.gz
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R
CHR="$SGE_TASK_ID"
OUT_BED_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Formatted/WM_M_musculus_sp_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed.gz

### RUN

Rscript --vanilla $R_SCRIPT \
$CASTANEUS_FILE \
$DOMESTICUS_FILE \
$MUSCULUS_FILE \
$FUNCTIONS_FILE \
$CHR \
$OUT_BED_FILE 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0