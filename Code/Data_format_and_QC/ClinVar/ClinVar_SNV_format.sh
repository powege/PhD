#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N ClinVar_format
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
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/ClinVar/ClinVar_SNV_format.R  \
"/well/lindgren/George/Data/ClinVar/raw/variant_summary.txt" \
"/well/lindgren/George/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
