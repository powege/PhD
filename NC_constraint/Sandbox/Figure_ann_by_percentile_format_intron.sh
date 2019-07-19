#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 16
#$ -N Intron_by_percentile
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

# human
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_annotation_by_percentile_format_intron_2.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/human_annotation_territory_by_percentile_650_50.csv" \
"human"

# mouse
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_annotation_by_percentile_format_intron_2.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/mouse_annotation_territory_by_percentile_650_50.csv" \
"mouse"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

