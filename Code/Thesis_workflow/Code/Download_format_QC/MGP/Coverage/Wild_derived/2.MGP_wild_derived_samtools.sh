#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-7 -tc 7
#$ -N MGP_wild_samtools
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
# parralel by number of strains in wildSTRAINS.txt
# run samtools for bam file
# splits into chr specific files with cols POS; DP
# gzip output: [strain]_samtools_DP_chr[chrom].txt.gz


### SET VARS

# set working directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Raw/
# set your samtools variable
SAM=/gpfs0/apps/well/samtools/1.4.1/bin/samtools
# set $SGE to file name
STRAIN_PREFIX=`sed -n -e "$SGE_TASK_ID p" "$PED_ROOT"wildSTRAINS.txt`


### SAMTOOLS

# run samtools to get sequence depth (-aa argument outputs all bp)
$SAM depth -aa "$PED_ROOT""$STRAIN_PREFIX"".bam" > "$PED_ROOT""$STRAIN_PREFIX".tmp && \
mv "$PED_ROOT""$STRAIN_PREFIX".tmp "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt


### FORMAT OUTPUT

# subset autosomes and X
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
$1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
$1 == "17" || $1 == "18" || $1 == "19" || $1 == "X")' "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt > "$PED_ROOT""$STRAIN_PREFIX".tmp && \
mv "$PED_ROOT""$STRAIN_PREFIX".tmp "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt

# split by chr
cd $PED_ROOT
awk '{print $0>$1 "_'$STRAIN_PREFIX'_samtools_DP.txt"}' "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt
for CHR in {1..19} X
do
  mv "$PED_ROOT""$CHR"_"$STRAIN_PREFIX"_samtools_DP.txt "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP_chr"$CHR".txt
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
rm "$PED_ROOT""$STRAIN_PREFIX".bam
rm "$PED_ROOT""$STRAIN_PREFIX"_samtools_DP.txt


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





