#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
# #$ -pe shmem 7
#$ -t 1-23 -tc 16
#$ -N alignment_format
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

### RUN

Rscript --vanilla /well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Ensembl/Alignment/2.pairwise_alignment_maf_to_bed_format_QC.R \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Raw/ \
/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/ \
hsap_grch38.v.mmus_grcm38.lastz_net/ \
hsap_grch38.v.mmus_grcm38.lastz_net. \
hsap_grch38_v_mmus_grcm38_v101_alignment_chr \
$SGE_TASK_ID


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
