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

# ARGS
# species <- args[1]
# chromo <- as.integer(args[2])
# window.size <- as.integer(args[3])
# window.shift <- as.integer(args[4])
# ref.file <- args[5]
# vcf.file <- args[6]
# mask.file <- args[7]
# alignment.file <- args[8]
# mu_table.file <- args[9]
# out.file <- args[10]

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Constraint_by_window/Constraint_variables_by_window.R \
mouse \
"$SGE_TASK_ID" \
750 \
50 \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/vcf_QCed_VEP/MGP_QCed_VEP_v94_hom_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/Masks/MGP_Mask_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/Ensembl/Alignment/Formatted/M_HtoM_alignment_short.txt \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table \
/well/lindgren/George/Data/NC_constraint/Constraint/mouse_constraint_variables_by_window_chr"$SGE_TASK_ID"_750_50.csv

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

