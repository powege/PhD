#!/bin/bash

#$ -N k7mer_SNV_rates
#$ -P lindgren.prjc -q short.qc
#$ -cwd -V
#$ -pe shmem 5
#$ -t 1-22 -tc 5
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

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_counts_by_chr.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/vcf_QCed_VEP/MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_counts_unmasked_chr"$SGE_TASK_ID".table \
"$SGE_TASK_ID" \
unmasked

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_counts_by_chr.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/MGP/vcf_QCed_VEP/MGP_v5_snps_QCed_hom_chr"$SGE_TASK_ID".vcf \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_counts_chr"$SGE_TASK_ID".table \
"$SGE_TASK_ID" \
all

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_counts_by_chr.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/1KGP/vcf_QCed_VEP/1KGP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/NC_constraint/SNV_rates/1KGP_7mer_SNV_counts_unmasked_chr"$SGE_TASK_ID".table \
"$SGE_TASK_ID" \
unmasked

Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_counts_by_chr.R \
/well/lindgren/George/Data/Ensembl/Reference/Formatted/Human_REF_sm_Ensembl_GRCm38_v94_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/1KGP/vcf_QCed_VEP/1KGP_phase3_QCed_VEP_v94_allPASS_chr"$SGE_TASK_ID".txt \
/well/lindgren/George/Data/NC_constraint/SNV_rates/1KGP_7mer_SNV_counts_chr"$SGE_TASK_ID".table \
"$SGE_TASK_ID" \
all

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
