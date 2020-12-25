#!/bin/bash

########
### SET VARS
########

# set working directory for files to download into
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Raw/

########
### download softmasked refernece fasta by chromosome
########

wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.chromosome.*.fa.gz \
-P "$PED_ROOT"

wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.*.fa.gz \
-P "$PED_ROOT"