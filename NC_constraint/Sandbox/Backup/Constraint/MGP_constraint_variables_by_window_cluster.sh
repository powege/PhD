#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-19 -tc 16
#$ -N M_constraint_vars
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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Constraint/MGP_constraint_variables_by_window_cluster.R "$SGE_TASK_ID"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

