#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-22 -tc 16
#$ -N WM_musculus_samtools
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

### SCRIPT that:

### SET VARS

# set working directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/
# set your samtools variable
SAM=/gpfs0/apps/well/samtools/1.4.1/bin/samtools
# set $SGE to file name
STRAIN_PREFIX=`sed -n -e "$SGE_TASK_ID p" "$PED_ROOT"musculus_bam.txt`


### SAMTOOLS

# run samtools to get sequence depth (-aa argument outputs all bp)
$SAM depth -aa "$PED_ROOT""$STRAIN_PREFIX" > "$PED_ROOT""$STRAIN_PREFIX".tmp && \
mv "$PED_ROOT""$STRAIN_PREFIX".tmp "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt


### FORMAT OUTPUT

# subset autosomes and X
awk '($1 == "chr1"  || $1 == "chr2" || $1 == "chr3"  || $1 == "chr4" || $1 == "chr5" || $1 == "chr6"  || $1 == "chr7" || $1 == "chr8" ||\
$1 == "chr9"  || $1 == "chr10" || $1 == "chr11" || $1 == "chr12"  || $1 == "chr13" || $1 == "chr14" || $1 == "chr15" || $1 == "chr16"  ||\
$1 == "chr17" || $1 == "chr18" || $1 == "chr19" || $1 == "chrX")' "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt > "$PED_ROOT""$STRAIN_PREFIX".tmp && \
mv "$PED_ROOT""$STRAIN_PREFIX".tmp "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt

# split by chr
cd $PED_ROOT
awk '{print $0>$1 "_'$STRAIN_PREFIX'_samtools_DP.txt"}' "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt
for CHR in {1..19} X
do
  mv "$PED_ROOT"chr"$CHR"_"$STRAIN_PREFIX"_samtools_DP.txt "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP_chr"$CHR".txt
done

# subset POS and DP cols
for CHR in {1..19} X
do
  awk '{print $2, $3}' "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP_chr"$CHR".txt > "$PED_ROOT""$STRAIN_PREFIX""$CHR".tmp && \
  mv "$PED_ROOT""$STRAIN_PREFIX""$CHR".tmp "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP_chr"$CHR".txt
done

# gzip
gzip "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP_chr*.txt

# remove intermediate files
# rm "$PED_ROOT""$STRAIN_PREFIX"
# rm "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"

exit 0


# # check file lengths
# for FILE_NO in {1..35}
# do
#   # set $infile1 to file name
#   infile1=`sed -n -e "$FILE_NO p" "$PED_ROOT"FILES`
#   # check file length
#   wc -l "$PED_ROOT""$infile1"_depth.txt 
# done





