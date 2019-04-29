#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -N Annotation_format
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
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Ensembl/Ensembl_annotation_format.R  \
"/well/lindgren/George/Data/Ensembl/Annotation/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff" \
"/well/lindgren/George/Data/Ensembl/Annotation/Homo_sapiens.GRCh38.94.gtf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv"

# run Ensembl_annotation_format.R for mouse
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Ensembl/Ensembl_annotation_format.R  \
"/well/lindgren/George/Data/Ensembl/Annotation/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mus_musculus.GRCm38.94.gtf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv"

# run Ensembl_unannotated_format.R for human
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Ensembl/Ensembl_unannotated_format.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_unannotated.csv" \
homo_sapien

# run Ensembl_unannotated_format.R for mouse
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Ensembl/Ensembl_unannotated_format.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_unannotated.csv" \
mus_musculus

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
