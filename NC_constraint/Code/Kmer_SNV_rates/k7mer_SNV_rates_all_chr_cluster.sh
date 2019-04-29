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

# Mouse all
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_rates_all_chr.R \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_counts_chr \
19 \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table

# Mouse unmasked
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_rates_all_chr.R \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_counts_unmasked_chr \
19 \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates_unmasked.table

# Human all
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_rates_all_chr.R \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_counts_chr \
22 \
/well/lindgren/George/Data/NC_constraint/SNV_rates/MGP_7mer_SNV_rates.table

# Mouse unmasked
Rscript --vanilla /well/lindgren/George/Workflows/NC_constraint/Code/Kmer_SNV_rates/k7mer_SNV_rates_all_chr.R \
/well/lindgren/George/Data/NC_constraint/SNV_rates/1KGP_7mer_SNV_counts_unmasked_chr \
22 \
/well/lindgren/George/Data/NC_constraint/SNV_rates/1KGP_7mer_SNV_rates_unmasked.table

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
