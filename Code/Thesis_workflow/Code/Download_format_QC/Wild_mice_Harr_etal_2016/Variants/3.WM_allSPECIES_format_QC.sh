#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -N WM_allSPECIES_vcf_QC
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

# directories
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/
FORMATTED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Formatted/
# files
IN_FILE=AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf
OUT_PREFIX=WM_Harr_etal_2016_allSPECIES_snps_PASS_chr
# R script
R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Wild_mice_Harr_etal_2016/Variants/2.WM_vcf_format_QC.R

### BASH FORMATTING

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE".gz > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE"

### RUN R SCRIPT

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE" \
"$FORMATTED_ROOT""$OUT_PREFIX"

# remove tmp files
rm "$RAW_ROOT""$IN_FILE"

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0
