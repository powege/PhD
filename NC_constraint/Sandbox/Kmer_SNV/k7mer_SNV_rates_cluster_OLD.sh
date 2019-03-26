#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -t 1-19 -tc 16
#$ -N k7mer_SNV_rates
#$ -o /well/lindgren/George/Workflows/NC_constraint/Log/
#$ -e /well/lindgren/George/Workflows/NC_constraint/Log/

#############
### Error: cannot allocate vector of size 1.5 Gb
#############

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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV/k7mer_SNV_rates_cluster.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/vcf_QCed_VEP/MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_chr"$SGE_TASK_ID".table

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
