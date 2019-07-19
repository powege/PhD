#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N Alignment_MAF_format
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

# run Ensembl_alignment_MAF_format.R
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Alignment/Ensembl_alignment_MAF_format.R  \
"/well/lindgren/George/Data/Ensembl/Alignment/Raw/" \
"/well/lindgren/George/Data/Ensembl/Alignment/Formatted/" \
"homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/" \
"homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr" \
"H_HtoM_alignment_long_chr" 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
