#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-22 -tc 16
#$ -N StrictMask_fasta_format
#$ -o /well/lindgren/George/Workflows/NC_constraint/Log/
#$ -e /well/lindgren/George/Workflows/NC_constraint/Log/

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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Masking/1KGP_mask_format.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/1KGP/StrictMask/Raw/20160622.chr"$SGE_TASK_ID".mask.fasta \
/well/lindgren/George/Data/1KGP/StrictMask/Formatted/1KGP_Mask_chr"$SGE_TASK_ID".txt 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
