#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -N methylation_format
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
R_SCRIPT=/well/lindgren/George/Code/Workflows/Sand/Ensembl_methylation_format_QC.R

CG_FILE=
ES_FILE=
NPC_FILE=
OUT_FILE=

### RUN

Rscript --vanilla $R_SCRIPT \
$CG_FILE \
$ES_FILE \
$NPC_FILE \
$OUT_FILE

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0