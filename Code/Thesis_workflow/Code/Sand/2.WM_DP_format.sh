#!/bin/bash

#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -N Wild_coverage_format1
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

# Script that concatinates coverage files by chr

# set working directory
PED_ROOT=/well/lindgren/George/Data/Wild_mice/coverage/raw/

# concatenate files
cd $PED_ROOT
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
cat "$PED_ROOT"*_depth_less10X_chr"$CHR".txt > "$PED_ROOT"Wild_mice_depth_less10X_chr"$CHR".txt
echo "$CHR"
done