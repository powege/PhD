#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 7
#$ -N Ann_by_per
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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Figures/Figure_annotations_by_percentile.R

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
