#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -t 1-19 -tc 16
#$ -N MGP_coverage_format2
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

# Script that pastes RD for eaach strain to one file

# set working directory
PED_ROOT=/well/lindgren/George/Data/MGP/bam/

# subset CHR POS form file as template
awk '{print $1, $2}' "$PED_ROOT"AKR_J.bam_"$SGE_TASK_ID".txt > "$PED_ROOT"MGP_coverage_allPOS_chr"$SGE_TASK_ID".txt

for FILE_NO in {1..40}
do

  # set $infile1 to file name
  infile1=`sed -n -e "$FILE_NO p" "$PED_ROOT"FILES`

  # check file length
  # wc -l "$PED_ROOT""$infile1".txt

  # subset read depth
  awk '{print $3}' "$PED_ROOT""$infile1"_"$SGE_TASK_ID".txt > "$PED_ROOT"tmp_RD_"$SGE_TASK_ID"

  # paste to file
  paste "$PED_ROOT"MGP_coverage_allPOS_chr"$SGE_TASK_ID".txt "$PED_ROOT"tmp_RD_"$SGE_TASK_ID" > "$PED_ROOT"tmp_file_"$SGE_TASK_ID"
  mv "$PED_ROOT"tmp_file_"$SGE_TASK_ID" "$PED_ROOT"MGP_coverage_allPOS_chr"$SGE_TASK_ID".txt

  echo $FILE_NO
  
done





#####

# for FILE_NO in {1..40}
# do
# 
#   # set $infile1 to file name
#   infile1=`sed -n -e "$FILE_NO p" "$PED_ROOT"FILES`
#   # check file length
#   wc -l "$PED_ROOT""$infile1".txt
# 
# done


