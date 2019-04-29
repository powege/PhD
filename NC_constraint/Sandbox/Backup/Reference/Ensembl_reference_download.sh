#!/bin/bash

### Script that downloads and unzips Mouse reference genome (release 94) from Ensembl ftp

# Set working directory for files to download into
PED_ROOT=/well/lindgren/George/Data/Ensembl/Reference/Raw/

# Download mouse refernece fasta by chromosome
wget ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.*.fa.gz \
-P "$PED_ROOT"
# remove unwanted files
rm Mus_musculus.GRCm38.dna_sm.chromosome.Y.fa.gz
rm Mus_musculus.GRCm38.dna_sm.chromosome.MT.fa.gz
# unzip mouse files
gunzip Mus_musculus.GRCm38.dna_sm.chromosome.*.fa.gz

# Download human refernece fasta by chromosome
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.*.fa.gz \
-P "$PED_ROOT"
# remove unwanted files
rm Homo_sapiens.GRCh38.dna_sm.chromosome.Y.fa.gz
rm Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa.gz
# unzip mouse files
gunzip Homo_sapiens.GRCh38.dna_sm.chromosome.*.fa.gz