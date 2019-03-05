#!/bin/bash

### Script that downloads and unzips Mouse reference genome (release 94) from Ensembl ftp

# Set working directory for files to download into
cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/Ensembl_reference/

# Download Mouse refernece fasta by chromosome
wget ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.*.fa.gz

# remove unwanted files
rm Mus_musculus.GRCm38.dna.chromosome.Y.fa.gz
rm Mus_musculus.GRCm38.dna.chromosome.MT.fa.gz

# unzip mouse files
gunzip Mus_musculus.GRCm38.dna.chromosome.*.fa.gz


