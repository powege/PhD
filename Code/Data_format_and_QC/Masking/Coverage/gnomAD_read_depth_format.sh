#!/bin/bash
#$ -cwd -V
#$ -P lindgren.prjc -q short.qc
#$ -N gnomAD_coverage_format
#$ -o /well/lindgren/George/Log/
#$ -e /well/lindgren/George/Log/

### Script that formats gnomAD coverage file to output CHR and POS with < 90% of individuals have min 10X coverage:

# Set variables
PED_ROOT=/well/lindgren/George/Data/gnomAD/coverage/
IN_FILE=release_2.1_coverage_genomes_gnomad.genomes.coverage.summary.tsv
OUT_FILE=gnomad_v2.1_fraction90_depth_less10X_chr"$CHR".tsv 

# unzip 
gunzip "$PED_ROOT""$IN_FILE".gz
echo "file unzipped"

# subset autosome POS that have < 90% of individauls with min 10X coverage
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
 $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
  $1 == "17" || $1 == "18" || $1 == "19" || $1 == "20" || $1 == "21"|| $1 == "22")\
   && ($7<0.9)' "$PED_ROOT""$IN_FILE" > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT"tmp_file
   
# subset CHR and POS
awk '{print $1, $2}' "$PED_ROOT"tmp_file > "$PED_ROOT"tmp
mv "$PED_ROOT"tmp "$PED_ROOT"tmp_file

# split file by chr
cd $PED_ROOT
awk '{print $0>$1 "_tmp_file.tsv"}' "$PED_ROOT"tmp_file
for CHR in {1..22}
do
mv "$CHR"_tmp_file.tsv "$OUT_FILE""$CHR".tsv
done

# subset POS
for CHR in {1..22}
do
awk '{print $2}' "$PED_ROOT""$OUT_FILE""$CHR".tsv > "$PED_ROOT"tmp"$CHR".tsv
mv "$PED_ROOT"tmp"$CHR".tsv "$PED_ROOT""$OUT_FILE""$CHR".tsv
done

# zip 
gzip "$PED_ROOT""$IN_FILE"
echo "file zipped"










   