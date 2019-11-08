#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N Annotation_overlap_summary
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

# HUMAN
Rscript --vanilla /well/lindgren/George/Code/Interspecific_mapping/Annotations/HM_annotation_overlap_summary.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr" \
"/well/lindgren/George/Data/Ensembl/Alignment/Formatted/H_annotation_overlap_summary.csv" \
"human"

# MOUSE
Rscript --vanilla /well/lindgren/George/Code/Interspecific_mapping/Annotations/HM_annotation_overlap_summary.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_chr" \
"/well/lindgren/George/Data/Ensembl/Alignment/Formatted/M_annotation_overlap_summary.csv" \
"mouse"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
