#!/bin/bash

# Script that:
# Downloads .bam file from ftp for each strain
# Run samtools to get sequence depth for all bp
# Subsets CHR and POS for all bases with coverage < 10x in autosomes and X chromosome

# Set working directory
PED_ROOT=/NGS/working_projects/MGP_SD/
SAM=/NGS/Software/samtools-1.3.1/samtools

# get a list of .bam files in ftp
curl -l ftp://ftp-mouse.sanger.ac.uk/current_bams/ > "$PED_ROOT"FILES
grep ".bam$" "$PED_ROOT"FILES > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"FILES


for bam in `cat "$PED_ROOT"FILES`
do

# Download gnomAD genomes vcf from ftp
wget -P "$PED_ROOT" ftp://ftp-mouse.sanger.ac.uk/current_bams/$bam 

# run samtools to get sequence depth (-aa argument outputs all bp)
$SAM depth -aa "$PED_ROOT"$bam > "$PED_ROOT"$bam.txt

# Subset all positions on autosomes and X chromosome with coverage < 10x
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
 $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
  $1 == "17" || $1 == "18" || $1 == "19"  || $1 == "X")\
   && ($3<10)' "$PED_ROOT"$bam.txt > "$PED_ROOT"$bam.tmp && mv "$PED_ROOT"$bam.tmp "$PED_ROOT"$bam.txt

# subset CHR and POS
awk '{print $1, $2}' "$PED_ROOT"$bam.txt > "$PED_ROOT"$bam.tmp && mv "$PED_ROOT"$bam.tmp "$PED_ROOT"$bam.txt

# remove .bam file
rm "$PED_ROOT"$bam

done

# concatenate all files and remove duplicates 

# delete raw files



