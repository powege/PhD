#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-22 -tc 16
#$ -N annnotation_model_format
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
R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Ensembl/Annotation/2.annotation_model_format_QC.R 

MOUSE_REG_HEART=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/mus_musculus.GRCm38.heart_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz
MOUSE_REG_KIDNEY=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/mus_musculus.GRCm38.kidney_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz
MOUSE_REG_SPLEEN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/mus_musculus.GRCm38.spleen_adult_8_weeks.Regulatory_Build.regulatory_activity.20180516.gff.gz
MOUSE_GENE_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/Mus_musculus.GRCm38.101.gtf.gz
MOUSE_N_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz
MOUSE_SM_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz
MOUSE_BP_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_summary.bed
MOUSE_TAD_IN=/well/lindgren/George/Data/Thesis_workflow/Data/TADs/Formatted/Mouse_TAD_boundries_liver_ESC_25kb_GRCm38.csv 
MOUSE_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz

HUMAN_REG_HEART=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/homo_sapiens.GRCh38.heart.Regulatory_Build.regulatory_activity.20190329.gff.gz
HUMAN_REG_KIDNEY=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/homo_sapiens.GRCh38.kidney.Regulatory_Build.regulatory_activity.20190329.gff.gz
HUMAN_REG_SPLEEN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/homo_sapiens.GRCh38.spleen.Regulatory_Build.regulatory_activity.20190329.gff.gz
HUMAN_GENE_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Raw/Homo_sapiens.GRCh38.101.gtf.gz
HUMAN_N_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_N.bed.gz
HUMAN_SM_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_sm.bed.gz
HUMAN_BP_IN=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_summary.bed
HUMAN_TAD_IN=/well/lindgren/George/Data/Thesis_workflow/Data/TADs/Formatted/Human_TAD_boundries_liver_ESC_25kb_GRCh38.csv 
HUMAN_OUT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz

### RUN

# Mouse 
Rscript --vanilla $R_SCRIPT \
$MOUSE_REG_HEART \
$MOUSE_REG_KIDNEY \
$MOUSE_REG_SPLEEN \
$MOUSE_GENE_IN \
$MOUSE_TAD_IN \
$MOUSE_BP_IN \
$MOUSE_N_IN \
$MOUSE_SM_IN \
$MOUSE_OUT \
$SGE_TASK_ID

# Human
Rscript --vanilla $R_SCRIPT \
$HUMAN_REG_HEART \
$HUMAN_REG_KIDNEY \
$HUMAN_REG_SPLEEN \
$HUMAN_GENE_IN \
$HUMAN_TAD_IN \
$HUMAN_BP_IN \
$HUMAN_N_IN \
$HUMAN_SM_IN \
$HUMAN_OUT \
$SGE_TASK_ID


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0