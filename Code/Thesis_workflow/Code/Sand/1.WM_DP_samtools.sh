#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-24 -tc 16
#$ -N Wild_samtools
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

# Script that:
# run samtools to get sequence depth for all bp
# subsets autosomes
# subsets CHR POS with read depth < 10 
# splits file by chr

# set working directory
PED_ROOT=/well/lindgren/George/Data/Wild_mice/coverage/raw/

#############
### SAMTOOLS
###########

#  make sure you have the right $MODULEPATH so the command can find everything
# module use -a /mgmt/modules/eb/modules/all

# load the samtools module
# module load SAMtools/1.8-intel-2018a

# set your samtools variable
# SAM=$EBROOTSAMTOOLS/bin/samtools
SAM=/gpfs0/apps/well/samtools/1.4.1/bin/samtools

# set $SGE to file name
infile1=`sed -n -e "$SGE_TASK_ID p" "$PED_ROOT"FILES`
# infile1=14_sorted.cigar.nodups.realigned.recalibrated.sorted.bam

# run samtools to get sequence depth (-aa argument outputs all bp)
$SAM depth -aa "$PED_ROOT""$infile1" > "$PED_ROOT""$infile1"_depth.txt 


###########
### FORMAT
#########

# subset autosomes
awk '($1 == "chr1"  || $1 == "chr2" || $1 == "chr3"  || $1 == "chr4" || $1 == "chr5" || $1 == "chr6"  || $1 == "chr7" || $1 == "chr8" ||\
$1 == "chr9"  || $1 == "chr10" || $1 == "chr11" || $1 == "chr12"  || $1 == "chr13" || $1 == "chr14" || $1 == "chr15" || $1 == "chr16"  ||\
$1 == "chr17" || $1 == "chr18" || $1 == "chr19")' "$PED_ROOT""$infile1"_depth.txt > "$PED_ROOT""$infile1"tmp_file
mv "$PED_ROOT""$infile1"tmp_file "$PED_ROOT""$infile1"_depth.txt 

# Subset all positions with coverage < 10x
awk '($3<10)' "$PED_ROOT""$infile1"_depth.txt > "$PED_ROOT""$infile1"_depth_less10X.txt

# subset CHR and POS
awk '{print $1, $2}' "$PED_ROOT""$infile1"_depth_less10X.txt > "$PED_ROOT""$infile1".tmp
mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1"_depth_less10X.txt

# split by chr
cd $PED_ROOT
awk '{print $0>$1 "_'$infile1'_depth_less10X.txt"}' "$PED_ROOT""$infile1"_depth_less10X.txt
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
  mv chr"$CHR"_"$infile1"_depth_less10X.txt "$infile1"_depth_less10X_chr"$CHR".txt
done

# remove 
rm "$PED_ROOT""$infile1"_depth.txt 
rm "$PED_ROOT""$infile1"_depth_less10X.txt




##########
# # Subsets CHR and POS for all bases with coverage < 10x in autosomes

# # Subset all positions on autosomes with coverage < 10x
# awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
#  $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
#   $1 == "17" || $1 == "18" || $1 == "19")\
#    && ($3<10)' "$PED_ROOT""$infile1".txt > "$PED_ROOT""$infile1".tmp && mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1".txt
   
# # Subset all positions on autosomes
# awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
#  $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
#   $1 == "17" || $1 == "18" || $1 == "19")' "$PED_ROOT""$infile1".txt > "$PED_ROOT""$infile1".tmp && mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1".txt 
   
# # subset CHR and POS
# awk '{print $1, $2}' "$PED_ROOT""$infile1".txt > "$PED_ROOT""$infile1".tmp && mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1".txt 

# delete raw files


# # check file lengths
# for FILE_NO in {1..35}
# do
#   # set $infile1 to file name
#   infile1=`sed -n -e "$FILE_NO p" "$PED_ROOT"FILES`
#   # check file length
#   wc -l "$PED_ROOT""$infile1"_depth.txt 
# done





