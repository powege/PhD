#!/bin/bash

#$ -cwd -V
#$ -l h_rt=24:00:00
#$ -P lindgren.prjb -q long.qb
#$ -t 1-40 -tc 16
#$ -N MGP_bam


# Script that:
# Run samtools to get sequence depth for all bp
# Subsets CHR and POS for all bases with coverage < 10x in autosomes and X chromosome

# set working directory
PED_ROOT=/well/lindgren/George/Workflows/NC_constraint/Data/MGP/bam/

#  make sure you have the right $MODULEPATH so the command can find everything
module use -a /mgmt/modules/eb/modules/all

# load the samtools module
module load SAMtools/1.8-intel-2018a

# set your samtools variable
SAM=$EBROOTSAMTOOLS/bin/samtools

# set $SGE to file name
infile1=`sed -n -e "$SGE_TASK_ID p" "$PED_ROOT"FILES`

# run samtools to get sequence depth (-aa argument outputs all bp)
$SAM depth -aa "$PED_ROOT""$infile1" > "$PED_ROOT""$infile1".txt 

# Subset all positions on autosomes and X chromosome with coverage < 10x
awk '($1 == "1"  || $1 == "2" || $1 == "3"  || $1 == "4" || $1 == "5" || $1 == "6"  || $1 == "7" || $1 == "8" ||\
 $1 == "9"  || $1 == "10" || $1 == "11" || $1 == "12"  || $1 == "13" || $1 == "14" || $1 == "15" || $1 == "16"  ||\
  $1 == "17" || $1 == "18" || $1 == "19"  || $1 == "X")\
   && ($3<10)' "$PED_ROOT""$infile1".txt > "$PED_ROOT""$infile1".tmp && mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1".txt 
   
# subset CHR and POS
awk '{print $1, $2}' "$PED_ROOT""$infile1".txt > "$PED_ROOT""$infile1".tmp && mv "$PED_ROOT""$infile1".tmp "$PED_ROOT""$infile1".txt 

# delete raw files





