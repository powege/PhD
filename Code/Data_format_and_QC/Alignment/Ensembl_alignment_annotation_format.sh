#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-23 -tc 16
#$ -N Alignment_ann_MAF_format
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

# run Ensembl_alignment_MAF_format.R
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Alignment/Ensembl_alignment_annotation_format.R  \
"/well/lindgren/George/Data/Ensembl/Alignment/Raw/" \
"/well/lindgren/George/Data/Ensembl/Alignment/Formatted/" \
"homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net/" \
"homo_sapiens_GRCh38_vs_mus_musculus_GRCm38_lastz_net.chr" \
"HM_alignment_annotation_long_Hchr" \
"HM_alignment_annotation_short_Hchr" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv" \
"$SGE_TASK_ID"


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
