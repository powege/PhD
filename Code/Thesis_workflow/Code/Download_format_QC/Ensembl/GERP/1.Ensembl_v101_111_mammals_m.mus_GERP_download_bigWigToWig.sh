#!/bin/bash

### LOAD MODULES
module load Kent_tools/401-gompi-2019b

### SET VARS
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/GERP/Raw/ # set working directory for files to download into
FORMATTED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/GERP/Formatted/

### download methylation big bed files
wget ftp://ftp.ensembl.org/pub/release-101/compara/conservation_scores/111_mammals.gerp_conservation_score/gerp_conservation_scores.mus_musculus.GRCm38.bw \
-P "$RAW_ROOT"

### convert bigBed to bed
bigWigToWig "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38.bw "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38.wig

### remove # lines
grep -v '^#' "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38.wig > "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38.wig.nohash

### split by chromosome
cd $RAW_ROOT
awk '{print $0>$1 "_gerp_conservation_scores.mus_musculus.GRCm38.bed"}' "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38.wig.nohash
for CHR in {1..19} X
do
mv "$RAW_ROOT""$CHR"_gerp_conservation_scores.mus_musculus.GRCm38.bed "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38_chr"$CHR".bed
done

### zip
for CHR in {1..19} X
do
gzip "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38_chr"$CHR".bed
echo "$CHR"
done

### move files
mv "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38_chr*.bed.gz "$FORMATTED_ROOT"

### remove intermederies 
rm "$RAW_ROOT"*_gerp_conservation_scores.mus_musculus.GRCm38.bed
rm "$RAW_ROOT"gerp_conservation_scores.mus_musculus.GRCm38.*
