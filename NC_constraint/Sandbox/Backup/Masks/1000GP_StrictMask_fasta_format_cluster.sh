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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Masks/1000GP_StrictMask_fasta_format.R \
/well/lindgren/George/Data/1000GP/Masks/Raw/20160622.chr"$SGE_TASK_ID".mask.fasta \
/well/lindgren/George/Data/1000GP/Masks/Formatted/1000GP_StrictMask_chr"$SGE_TASK_ID".txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
