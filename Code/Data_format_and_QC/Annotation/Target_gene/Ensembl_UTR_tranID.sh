#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 16
#$ -t 1-22 -tc 16
#$ -N UTR_link
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

# human
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Annotation/Target_gene/Ensembl_UTR_tranID.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Homo_sapiens.GRCh38.94.gtf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Human_UTR_tranID.csv"

# mouse
Rscript --vanilla /well/lindgren/George/Code/Data_format_and_QC/Annotation/Target_gene/Ensembl_UTR_tranID.R \
"/well/lindgren/George/Data/Ensembl/Annotation/Mus_musculus.GRCm38.94.gtf" \
"/well/lindgren/George/Data/Ensembl/Annotation/Mouse_UTR_tranID.csv"


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
