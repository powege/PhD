#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N Annotation_format
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
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Annotation/Format/Ensembl_WG_annotation_format_all.R \
"/well/lindgren/George/Data/Ensembl/Annotation/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff" \
"/well/lindgren/George/Data/Ensembl/Annotation/Homo_sapiens.GRCh38.94.gtf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID" \
"homo_sapien"

# mouse
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Annotation/Format/Ensembl_WG_annotation_format_all.R \
"/well/lindgren/George/Data/Ensembl/Annotation/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mus_musculus.GRCm38.94.gtf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_chr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID" \
"mus_musculus"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
