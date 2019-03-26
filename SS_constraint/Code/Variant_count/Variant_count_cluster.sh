#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -t 1-100 -tc 16
#$ -N Variant_count
#$ -o /well/lindgren/George/Workflows/SS_constraint/Log/
#$ -e /well/lindgren/George/Workflows/SS_constraint/Log/

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

Rscript --vanilla /well/lindgren/George/Workflows/SS_constraint/Code/Variant_count/Variant_count.R "$SGE_TASK_ID"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

#cat SS_gnomAD_variant_counts_chr1.csv SS_gnomAD_variant_counts_chr2.csv \
#SS_gnomAD_variant_counts_chr3.csv SS_gnomAD_variant_counts_chr4.csv \
#SS_gnomAD_variant_counts_chr5.csv SS_gnomAD_variant_counts_chr6.csv \ 
#SS_gnomAD_variant_counts_chr7.csv SS_gnomAD_variant_counts_chr8.csv \
#SS_gnomAD_variant_counts_chr9.csv SS_gnomAD_variant_counts_chr10.csv \
#SS_gnomAD_variant_counts_chr11.csv SS_gnomAD_variant_counts_chr12.csv \
#SS_gnomAD_variant_counts_chr13.csv SS_gnomAD_variant_counts_chr14.csv \
#SS_gnomAD_variant_counts_chr15.csv SS_gnomAD_variant_counts_chr16.csv \
#SS_gnomAD_variant_counts_chr17.csv SS_gnomAD_variant_counts_chr18.csv \
#SS_gnomAD_variant_counts_chr19.csv SS_gnomAD_variant_counts_chr20.csv \
#SS_gnomAD_variant_counts_chr21.csv SS_gnomAD_variant_counts_chr22.csv > SS_gnomAD_variant_counts.csv
