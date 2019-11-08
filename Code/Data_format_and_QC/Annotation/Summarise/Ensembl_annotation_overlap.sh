#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N annotation_overlap_summary
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

# Human
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Annotation/Ensembl_annotation_overlap_summary.R \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr""$SGE_TASK_ID"".txt" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_annotation_overlap_chr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID"

# Mouse
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Annotation/Ensembl_annotation_overlap_summary.R \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr""$SGE_TASK_ID"".txt" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_annotation_overlap_chr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
