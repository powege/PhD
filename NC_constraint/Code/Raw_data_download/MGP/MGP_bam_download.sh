#!/bin/bash

# Script that downloads .bam file from ftp for each strain

# Set working directory
PED_ROOT=/well/lindgren/George/Workflows/NC_constraint/Data/MGP/bam/


# get a list of .bam files in ftp
curl -l ftp://ftp-mouse.sanger.ac.uk/current_bams/ > "$PED_ROOT"FILES
grep ".bam$" "$PED_ROOT"FILES > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"FILES

for bam in `cat "$PED_ROOT"FILES`
do
# Download .bam files from ftp
wget -P "$PED_ROOT" ftp://ftp-mouse.sanger.ac.uk/current_bams/"$bam" 
done
