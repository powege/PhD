#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 3
#$ -t 1-22 -tc 16
#$ -N REF_fasta_format
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

module load R/3.4.3

R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Ensembl/Reference/2.GRC38_sm_reference_format.R
CHR="$SGE_TASK_ID"
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R

M_IN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Raw/Mus_musculus.GRCm38.dna_sm.chromosome."$SGE_TASK_ID".fa.gz
M_OUT_REF=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_seq_chr"$SGE_TASK_ID".txt.gz
M_OUT_SM=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz
M_OUT_N=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz
M_OUT_BP=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_summary.bed

H_IN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Raw/Homo_sapiens.GRCh38.dna_sm.chromosome."$SGE_TASK_ID".fa.gz
H_OUT_REF=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_seq_chr"$SGE_TASK_ID".txt.gz
H_OUT_SM=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_sm.bed.gz
H_OUT_N=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_N.bed.gz
H_OUT_BP=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_summary.bed

Rscript --vanilla $R_SCRIPT \
$M_IN_FILE \
$M_OUT_REF \
$M_OUT_SM \
$M_OUT_N \
$M_OUT_BP \
$CHR \
$FUNCTIONS_FILE 

Rscript --vanilla $R_SCRIPT \
$H_IN_FILE \
$H_OUT_REF \
$H_OUT_SM \
$H_OUT_N \
$H_OUT_BP \
$CHR \
$FUNCTIONS_FILE 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
