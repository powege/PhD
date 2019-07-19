#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-100 -tc 20
#$ -N Terriroty_by_percentile
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

# human common
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ann_by_percentile_format.R \
"$SGE_TASK_ID" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv" \
"/well/lindgren/George/Data/NC_constraint/Constraint/Constraint_by_window_human_950_50_MAF001.csv" \
"/well/lindgren/George/Data/NC_constraint/Figures_and_tables/human_ann_by_percentile_MAF001.csv" \
"human"

# human all
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ann_by_percentile_format.R \
"$SGE_TASK_ID" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv" \
"/well/lindgren/George/Data/NC_constraint/Constraint/Constraint_by_window_human_950_50.csv" \
"/well/lindgren/George/Data/NC_constraint/Figures_and_tables/human_ann_by_percentile.csv" \
"human"

# mouse hom
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Data_format_and_QC/Figures_and_Tables/Figure_ann_by_percentile_format.R \
"$SGE_TASK_ID" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation_ranked.csv" \
"/well/lindgren/George/Data/NC_constraint/Constraint/Constraint_by_window_mouse_950_50.csv" \
"/well/lindgren/George/Data/NC_constraint/Figures_and_tables/mouse_ann_by_percentile.csv" \
"mouse"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

