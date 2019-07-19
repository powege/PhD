#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-100 -tc 20
#$ -N POS_ID_percentile
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

# module use -a /apps/eb/skylake/modules/all
# module load R/3.5.1-foss-2018b-X11-20180604

module load R/3.4.3

# human
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Constraint_by_window/POS_ID_by_percentile.R \
"$SGE_TASK_ID" \
"human"

# mouse
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Constraint_by_window/POS_ID_by_percentile.R \
"$SGE_TASK_ID" \
"mouse"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

