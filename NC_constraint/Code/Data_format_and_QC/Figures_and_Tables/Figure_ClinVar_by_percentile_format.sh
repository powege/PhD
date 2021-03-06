#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N ClinVar_by_percentile
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

# run Ensembl_annotation_format.R for human
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ClinVar_by_percentile.R  \
"/well/lindgren/George/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf" \
"/well/lindgren/George/Data/ClinVar/formatted/ClinVar_human_to_mouse_alignment.txt" \
"/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/ClinVar_SNVs_by_percentile_650_50.csv" 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
