#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N ref_bp
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

# human
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Reference/Ensembl_reference_n_bp.R \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr" \
"/well/lindgren/George/Data/Ensembl/Reference/Human_REF_n_bp.csv" \
"/well/lindgren/George/Data/Ensembl/Reference/Human_REF_ATGC_POS.csv" \
"human"

# mouse
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Reference/Ensembl_reference_n_bp.R \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr" \
"/well/lindgren/George/Data/Ensembl/Reference/Mouse_REF_n_bp.csv" \
"/well/lindgren/George/Data/Ensembl/Reference/Mouse_REF_ATGC_POS.csv" \
"mouse"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
