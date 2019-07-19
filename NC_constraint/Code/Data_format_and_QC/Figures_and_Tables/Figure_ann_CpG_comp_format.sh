#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N annotation_CpG
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

# human
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ann_CpG_comp_format.R \
"$SGE_TASK_ID" \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv \
/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv \
/well/lindgren/George/Data/NC_constraint/Figures_and_tables/Figure_annotation_CpG_proportion_human.csv

# mouse
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ann_CpG_comp_format.R \
"$SGE_TASK_ID" \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv \
/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation_ranked.csv \
/well/lindgren/George/Data/NC_constraint/Figures_and_tables/Figure_annotation_CpG_proportion_mouse.csv

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

