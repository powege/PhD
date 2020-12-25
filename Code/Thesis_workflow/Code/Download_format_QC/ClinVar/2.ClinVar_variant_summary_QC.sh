#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -N ClinVar_QC
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

Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/ClinVar/2.ClinVar_variant_summary_QC.sh 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0