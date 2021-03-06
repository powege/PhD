#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N GWAS_DDC_alignment_null_sampling
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
Rscript --vanilla /well/lindgren/George/Code/Interspecific_mapping/SNVs/HM_SNV_alignment_null_sampling.R \
"/well/lindgren/George/Data/GWAS/formatted/GWAS_DDC_mouse_mapping_Hchr""$SGE_TASK_ID"".csv" \
"/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_long_Hchr""$SGE_TASK_ID"".csv" \
"/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr""$SGE_TASK_ID"".txt" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/GWAS/formatted/GWAS_DDC_mouse_mapping_null_sampling_Hchr""$SGE_TASK_ID"".csv" \
"$SGE_TASK_ID" \
"1000"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
