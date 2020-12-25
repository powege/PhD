#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 5
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
R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Ensembl/Methylation/2.Ensembl_v101_m.mus_methylation_format.R

CG_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRC38_dna_CG.bed.gz
ES_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Methylation/Raw/ES_5mC_Stadler2011_PMID22170606.bed.gz
NPC_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Methylation/Raw/NPC_5mC_Stadler2011_PMID22170606.bed.gz
OUT_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz

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