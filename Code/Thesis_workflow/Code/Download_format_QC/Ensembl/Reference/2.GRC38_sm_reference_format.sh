#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 4
#$ -t 1-22 -tc 16
#$ -N reference_fasta_format
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

R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Ensembl/Reference/2.GRC38_sm_reference_format.R
CHR="$SGE_TASK_ID"
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R

MOUSE_FASTA_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Raw/Mus_musculus.GRCm38.dna_sm.chromosome."$SGE_TASK_ID".fa.gz
MOUSE_N_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz
MOUSE_SM_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz
MOUSE_CG_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRC38_dna_CG.bed.gz
MOUSE_BP_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_summary.bed
MOUSE_REF_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr"$SGE_TASK_ID".txt.gz

HUMAN_FASTA_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Raw/Homo_sapiens.GRCh38.dna_sm.chromosome."$SGE_TASK_ID".fa.gz
HUMAN_N_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_N.bed.gz
HUMAN_SM_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_sm.bed.gz
HUMAN_CG_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_reference_Ensembl_GRC38_CG.bed.gz
HUMAN_BP_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_summary.bed
HUMAN_REF_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_seq_chr"$SGE_TASK_ID".txt.gz

### RUN

Rscript --vanilla $R_SCRIPT \
$MOUSE_FASTA_IN \
$CHR \
$FUNCTIONS_FILE \
$MOUSE_N_OUT \
$MOUSE_SM_OUT \
$MOUSE_CG_OUT \
$MOUSE_BP_OUT \
$MOUSE_REF_OUT

Rscript --vanilla $R_SCRIPT \
$HUMAN_FASTA_IN \
$CHR \
$FUNCTIONS_FILE \
$HUMAN_N_OUT \
$HUMAN_SM_OUT \
$HUMAN_CG_OUT \
$HUMAN_BP_OUT \
$HUMAN_REF_OUT

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0