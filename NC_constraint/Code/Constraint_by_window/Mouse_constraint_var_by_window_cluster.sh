#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-19 -tc 16
#$ -N Mouse_constraint_vars
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

# mouse
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Constraint_by_window/Constraint_variables_by_window.R \
"$SGE_TASK_ID" \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/vcf_QCed_VEP/MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/Masks/MGP_Mask_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table \
mouse

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

