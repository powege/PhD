#!/bin/bash

# Script that:
# Run samtools to get sequence depth for all bp
# Subsets CHR and POS for all bases with coverage < 10x in autosomes and X chromosome

# Set working directory
PED_ROOT=/NGS/working_projects/MGP_SD/
SAM=/NGS/Software/samtools-1.3.1/samtools


for bam in `cat "$PED_ROOT"FILES`
do

# run samtools to get sequence depth (-aa argument outputs all bp)
"$SAM" depth -aa "$bam" > "$bam".txt

# Subset all positions on autosomes and X chromosome with coverage < 10x
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
 $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
  $1 == "17" || $1 == "18" || $1 == "19"  || $1 == "X")\
   && ($3<10)' "$bam".txt > "$bam".tmp && mv "$bam".tmp "$bam".txt

# subset CHR and POS
awk '{print $1, $2}' "$bam".txt > "$bam".tmp && mv "$bam".tmp "$bam".txt

done

# concatenate all files and remove duplicates 

# delete raw files



