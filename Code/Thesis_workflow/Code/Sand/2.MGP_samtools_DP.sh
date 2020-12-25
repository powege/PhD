#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjb -q long.qb
#$ -t 1-36 -tc 16
#$ -N MGP_samtools
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

# Script that:
# run samtools to get sequence depth for all bp
# subsets autosomes and X
# subsets CHR POS with read depth < 10 
# splits file by chr
# outputs strain specific file

# set working directory
PED_ROOT=/well/lindgren/George/Data/MGP/bam/

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

# run samtools to get sequence depth (-aa argument outputs all bp)
$SAM depth -aa "$PED_ROOT""$infile1" > "$PED_ROOT""$infile1"_depth.txt 


###########
### FORMAT
#########

# subset autosomes and X
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
$1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
$1 == "17" || $1 == "18" || $1 == "19" || $1 == "X")' "$PED_ROOT""$infile1"_depth.txt > "$PED_ROOT""$infile1"tmp_file
mv "$PED_ROOT""$infile1"tmp_file "$PED_ROOT""$infile1"_depth.txt 

# Subset all positions with coverage < 10x
awk '($3<10)' "$PED_ROOT""$infile1"_depth.txt > "$PED_ROOT""$infile1"_depth_less10X.txt

# subset CHR and POS
awk '{print $1, $2}' "$PED_ROOT""$infile1"_depth_less10X.txt > "$PED_ROOT""$infile1".tmp && mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1"_depth_less10X.txt

# split by chr
cd $PED_ROOT
awk '{print $0>$1 "_'$infile1'_depth_less10X.txt"}' "$PED_ROOT""$infile1"_depth_less10X.txt
for CHR in {1..19} .. X
do
  mv "$CHR"_"$infile1"_depth_less10X.txt "$infile1"_depth_less10X_chr"$CHR".txt
done

# clear temporary files
rm "$PED_ROOT""$infile1"_depth.txt
rm "$PED_ROOT""$infile1"_depth_less10X.txt



##########


# # check file lengths
# for FILE_NO in {1..35}
# do
#   # set $infile1 to file name
#   infile1=`sed -n -e "$FILE_NO p" "$PED_ROOT"FILES`
#   # check file length
#   wc -l "$PED_ROOT""$infile1"_depth.txt 
# done





