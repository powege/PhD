#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 5
#$ -t 1-19 -tc 16
#$ -N WM_musculus_n_DP
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

R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Wild_mice_Harr_etal_2016/Coverage/3.DP_less10X.R

RAW_PATH=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/
PREFIX_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/musculus_bam.txt
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R
CHR="$SGE_TASK_ID"
OUT_LONG_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/WM_musculus_n_DP_less10X_chr"$SGE_TASK_ID".csv.gz
OUT_BED_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Formatted/WM_musculus_DP_less10X_more0.1_chr"$SGE_TASK_ID".bed.gz

### RUN

Rscript --vanilla $R_SCRIPT \
$RAW_PATH \
$PREFIX_FILE \
$FUNCTIONS_FILE \
$CHR \
$OUT_LONG_FILE \
$OUT_BED_FILE 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0