#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N Alignment_CinVar_format
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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ClinVar_alignment_format.R  \
"/well/lindgren/George/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf" \
"/well/lindgren/George/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_long_chr" \
"/well/lindgren/George/Data/ClinVar/formatted/ClinVar_human_to_mouse_alignment.txt" 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
