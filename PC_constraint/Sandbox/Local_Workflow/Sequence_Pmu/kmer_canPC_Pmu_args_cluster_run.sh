#!/bin/sh
#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 3
#$ -S /bin/bash
#$ -t 1-22
#$ -o /well/lindgren/George/Workflows/PC_constraint/Cluster/
#$ -e /well/lindgren/George/Workflows/PC_constraint/Cluster/

echo "Started at: "`date +"%m/%d/%Y %H:%M:%S"`

### set root path of R script
PED_ROOT=/well/lindgren/George/Workflows/PC_constraint/Paper/Code/Sequence_Pmu/

### set args
# 1. mutation probabiliy table file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table")
# 2. mutation coding table file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/AA_mutation_table.csv")
# 3. transcript sequence file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_canPC_SEQ_QCed.csv")
# 4. chromosome (eg 21)
# 5. output file path (eg "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k5mer_canPC_Pmu_chr21.csv")
arg1=/well/lindgren/George/Workflows/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table
arg2=/well/lindgren/George/Workflows/PC_constraint/Paper/Data/Mu_rates/AA_mutation_table.csv
arg3=/well/lindgren/George/Workflows/PC_constraint/Paper/Data/Ensembl/H_canPC_SEQ_QCed.csv
arg4="$SGE_TASK_ID"
arg5=/well/lindgren/George/Workflows/PC_constraint/Paper/Data/Sequence_Pmu/H_k5mer_canPC_Pmu_chr"$SGE_TASK_ID".csv

### run R script
Rscript --vanilla "PED_ROOT"k5mer_canPC_Pmu_byCHR.R arg1 arg2 arg3 arg4 arg5

echo "Ended at: "`date +"%m/%d/%Y %H:%M:%S"`
