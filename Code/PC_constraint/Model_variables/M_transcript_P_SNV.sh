#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -N M_transcript_P_SNV
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

# mouse
Rscript --vanilla /well/lindgren/George/Code/PC_constraint/Model_variables/Transcript_P_SNV.R \
"/well/lindgren/George/Data/Mu_rates/AA_mutation_table.csv" \
"/well/lindgren/George/Data/Ensembl/BioMart/QCed/Ensembl_v94_mouse_canPC_seq_QCpass.csv" \
"/well/lindgren/George/Data/PC_constraint/Ensembl_v94_mouse_canPC_QCpass_pMu.csv" \
"mouse"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
