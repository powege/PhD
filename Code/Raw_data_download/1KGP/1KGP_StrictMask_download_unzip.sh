#!/bin/bash

### Script that downloads and unzips 1KGP StrictMask

# set file directory
PED_ROOT=/well/lindgren/George/Data/1KGP/StrictMask/

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.chr"$CHR".mask.fasta.gz \
-P "$PED_ROOT"
gunzip "$PED_ROOT"20160622.chr"$CHR".mask.fasta.gz
mv "$PED_ROOT"20160622.chr"$CHR".mask.fasta "$PED_ROOT"1KGP_StrictMask_chr"$CHR".fasta
done