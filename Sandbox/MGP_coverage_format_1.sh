#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -t 1-35 -tc 16
#$ -N MGP_coverage_format1
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

# Script that subsets autosomes and splits files by CHR

# set working directory
PED_ROOT=/well/lindgren/George/Data/MGP/bam/

# set $SGE to file name
infile1=`sed -n -e "$SGE_TASK_ID p" "$PED_ROOT"FILES`

# subset autosomes
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
$1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
$1 == "17" || $1 == "18" || $1 == "19")' "$PED_ROOT""$infile1".txt > "$PED_ROOT""$infile1"tmp_file
mv "$PED_ROOT""$infile1"tmp_file "$PED_ROOT""$infile1".txt 

# split by chr
cd $PED_ROOT
awk '{print $0>$1 "_'$infile1'.txt"}' "$PED_ROOT""$infile1".txt 
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
  mv "$CHR"_"$infile1".txt "$infile1"_"$CHR".txt
done



# tail -99999999 BUB_BnJ.bam.txt > test.txt
# infile1=test



