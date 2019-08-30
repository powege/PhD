#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-19 -tc 16
#$ -N MGP_coverage_format2
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

Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Masking/MGP_read_depth_format_2.R \
/well/lindgren/George/Data/MGP/bam/MGP_depth_less10X_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/bam/MGP_fraction_depth_less10X_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/bam/MGP_fraction90_depth_less10X_chr"$SGE_TASK_ID".txt \
"$SGE_TASK_ID"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


