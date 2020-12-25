#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-22 -tc 16
#$ -N feature_alignnment_summary
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

### LOAD MODULE
module load R/3.4.3

### SET VARS 
R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Human_mouse_mapping/Analysis/1.feature_overlap_summary.R

HUMAN_ANN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_v101_whole_genome_features_multicell.csv.gz
MOUSE_ANN_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_v101_whole_genome_features_multicell.csv.gz
ALIGNMENT_FILE=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/hsap_grch38_v_mmus_grcm38_v101_alignment_chr"$SGE_TASK_ID".bed.gz
CHR="$SGE_TASK_ID"
FUNCTIONS_FILE=/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R
OUT_ALL_FILE=/well/lindgren/George/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_summary.csv
OUT_SPECIFIC_FILE=/well/lindgren/George/Data/Thesis_workflow/Results/Human_mouse_mapping/hsap_grch38_v_mmus_grcm38_v101_tissue_feature_alignment_summary.csv

### RUN
Rscript --vanilla $R_SCRIPT \
$HUMAN_ANN_FILE \
$MOUSE_ANN_FILE \
$ALIGNMENT_FILE \
$CHR \
$FUNCTIONS_FILE \
$OUT_ALL_FILE \
$OUT_SPECIFIC_FILE 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0