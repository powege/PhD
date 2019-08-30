#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-22 -tc 16
#$ -N Ensembl_canPC_fraction
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

# ### HUMAN gnomAD
# Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Transcript_QC/Ensembl_canPC_fraction_masked.R \
# /well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv \
# /well/lindgren/George/Data/gnomAD/coverage/gnomad_v2.1_fraction90_depth_less10X_chr"$SGE_TASK_ID".tsv \
# /well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_gnomad_mask_chr"$SGE_TASK_ID".csv \
# "$SGE_TASK_ID"

### HUMAN 1KGP
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Transcript_QC/Ensembl_canPC_fraction_masked.R \
/well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv \
/well/lindgren/George/Data/1KGP/StrictMask/Formatted/1KGP_StrictMask_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_1KGP_mask_chr"$SGE_TASK_ID".csv \
"$SGE_TASK_ID"

# ### MOUSE
# Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Transcript_QC/Ensembl_canPC_fraction_masked.R \
# /well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_pos.csv \
# /well/lindgren/George/Data/MGP/bam/MGP_fraction90_depth_less10X_chr"$SGE_TASK_ID".txt \
# /well/lindgren/George/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_MGP_mask_chr"$SGE_TASK_ID".csv \
# "$SGE_TASK_ID"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


