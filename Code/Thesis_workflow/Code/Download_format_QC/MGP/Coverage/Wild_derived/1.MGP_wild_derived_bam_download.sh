#!/bin/bash

# create files with wildSTRAINS 
# download bam files for all wildSTRAINS 

# Set working directory
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Raw/

# create file of wild-derived strains in vcf
echo 'CAST_EiJ
WSB_EiJ
PWK_PhJ
MOLF_EiJ
SPRET_EiJ
LEWES_EiJ
ZALENDE_EiJ' > "$RAW_ROOT"wildSTRAINS.txt

# download .bam files from ftp

# wild derived
for bam in `cat "$RAW_ROOT"wildSTRAINS.txt`
do
wget -P "$RAW_ROOT" ftp://ftp-mouse.sanger.ac.uk/current_bams/"$bam".bam
done

