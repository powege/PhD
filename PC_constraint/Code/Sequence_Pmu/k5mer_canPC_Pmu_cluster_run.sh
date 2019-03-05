#!/bin/sh
#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 3
#$ -S /bin/bash
#$ -o /well/lindgren/George/Workflows/PC_constraint/Cluster/
#$ -e /well/lindgren/George/Workflows/PC_constraint/Cluster/

echo "Started at: "`date +"%m/%d/%Y %H:%M:%S"`

### source and load R module
source /etc/profile.d/modules.sh
module load R/3.4.3

### run R script
Rscript --vanilla /well/lindgren/George/Workflows/PC_constraint/Code/Sequence_Pmu/k5mer_canPC_Pmu.R

echo "Ended at: "`date +"%m/%d/%Y %H:%M:%S"`

