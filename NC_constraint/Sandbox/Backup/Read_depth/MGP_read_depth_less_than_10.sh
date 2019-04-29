#!/bin/bash

# set working directory
cd /well/lindgren/George/Workflows/NC_constraint/Data/MGP/bam/
PED_ROOT=/well/lindgren/George/Workflows/NC_constraint/Data/MGP/bam/

# concatenate files
cat "$PED_ROOT"*.bam.txt > "$PED_ROOT"MGP_low_depth_POS.txt
#cat LP_J.bam.txt LG_J.bam.txt > MGP_low_depth_POS.txt

# remove duplicates
sort -u "$PED_ROOT"MGP_low_depth_POS.txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"MGP_low_depth_POS.txt

# split into files by CHR
awk '{print > "MGP_low_depth_POS_chr"$1".txt"}' "$PED_ROOT"MGP_low_depth_POS.txt

# delete raw files

