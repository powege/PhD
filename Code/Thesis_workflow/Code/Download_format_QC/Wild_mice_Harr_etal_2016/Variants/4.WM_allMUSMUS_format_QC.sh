#!/bin/bash

#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 5
#$ -N WM_allMUSMUS_vcf_QC
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
IN_FILE_Mmc_CAST=Mmc_CAST.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmd_FRA=Mmd_FRA.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmd_GER=Mmd_GER.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmd_HEL=Mmd_HEL.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmd_IRA=Mmd_IRA.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmm_AFG=Mmm_AFG.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmm_CZE=Mmm_CZE.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
IN_FILE_Mmm_KAZ=Mmm_KAZ.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf
OUT_PREFIX=WM_Harr_etal_2016_allMUSMUS_snps_PASS_chr

# R script
R_SCRIPT=/well/lindgren/George/Code/Thesis_workflow/Code/Download_format_QC/Wild_mice_Harr_etal_2016/Variants/2.WM_vcf_format_QC.R


### BASH FORMATTING

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmc_CAST".gz > "$RAW_ROOT""$IN_FILE_Mmc_CAST"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmc_CAST" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmc_CAST".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmc_CAST" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmc_CAST"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmc_CAST" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmc_CAST"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmd_FRA".gz > "$RAW_ROOT""$IN_FILE_Mmd_FRA"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmd_FRA" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_FRA".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmd_FRA" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_FRA"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmd_FRA" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_FRA"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmd_GER".gz > "$RAW_ROOT""$IN_FILE_Mmd_GER"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmd_GER" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_GER".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmd_GER" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_GER"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmd_GER" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_GER"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmd_HEL".gz > "$RAW_ROOT""$IN_FILE_Mmd_HEL"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmd_HEL" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_HEL".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmd_HEL" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_HEL"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmd_HEL" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_HEL"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmd_IRA".gz > "$RAW_ROOT""$IN_FILE_Mmd_IRA"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmd_IRA" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_IRA".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmd_IRA" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_IRA"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmd_IRA" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmd_IRA"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmm_AFG".gz > "$RAW_ROOT""$IN_FILE_Mmm_AFG"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmm_AFG" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_AFG".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmm_AFG" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_AFG"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmm_AFG" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_AFG"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmm_CZE".gz > "$RAW_ROOT""$IN_FILE_Mmm_CZE"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmm_CZE" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_CZE".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmm_CZE" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_CZE"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmm_CZE" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_CZE"

# unzip (keep original zipped file)
gunzip -c "$RAW_ROOT""$IN_FILE_Mmm_KAZ".gz > "$RAW_ROOT""$IN_FILE_Mmm_KAZ"
### remove comment lines (save as .meta file) This is necessary for R script
awk '/^#/' "$RAW_ROOT""$IN_FILE_Mmm_KAZ" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_KAZ".meta
awk '!/^#/' "$RAW_ROOT""$IN_FILE_Mmm_KAZ" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_KAZ"
### subset columns (CHROM POS ID REF ALT QUAL FILTER)
awk '{print $1, $2, $3, $4, $5, $6, $7}' "$RAW_ROOT""$IN_FILE_Mmm_KAZ" > "$RAW_ROOT"tmp && mv "$RAW_ROOT"tmp "$RAW_ROOT""$IN_FILE_Mmm_KAZ"


### RUN R SCRIPT

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmc_CAST" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmc_CAST_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmd_FRA" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmd_FRA_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmd_GER" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmd_GER_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmd_HEL" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmd_HEL_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmd_IRA" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmd_IRA_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmm_AFG" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmm_AFG_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmm_CZE" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmm_CZE_snps_PASS_pop_tmpfile_chr

Rscript --vanilla $R_SCRIPT \
"$RAW_ROOT""$IN_FILE_Mmm_KAZ" \
"$RAW_ROOT"WM_Harr_etal_2016_Mmm_KAZ_snps_PASS_pop_tmpfile_chr


### BASH FORMATTING

# unzip files
gunzip "$RAW_ROOT"*_pop_tmpfile_chr*.vcf.gz

# concatiante files
for CHR in {1..19} X
do
cat "$RAW_ROOT"*_pop_tmpfile_chr"$CHR".vcf > "$RAW_ROOT""$OUT_PREFIX""$CHR".vcf
done

# remove duplicates 
for CHR in {1..19} X
do
awk '!a[$0]++' "$RAW_ROOT""$OUT_PREFIX""$CHR".vcf > "$RAW_ROOT""$OUT_PREFIX""$CHR".vcf.tmp && mv "$RAW_ROOT""$OUT_PREFIX""$CHR".vcf.tmp "$RAW_ROOT""$OUT_PREFIX""$CHR".vcf
done

# zip
gzip "$RAW_ROOT""$OUT_PREFIX"*.vcf

# move dir
mv "$RAW_ROOT""$OUT_PREFIX"*.vcf.gz "$FORMATTED_ROOT"

# remove tmp files
rm "$RAW_ROOT"*.vcf


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0

