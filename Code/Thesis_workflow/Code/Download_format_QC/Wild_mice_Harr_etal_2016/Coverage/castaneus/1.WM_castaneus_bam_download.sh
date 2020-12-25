#!/bin/bash

# Set working directory
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/

# create file of  bam files

# Mus castaneus
# 10 individuals
echo 'H12_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H14_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H15_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H24_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H26_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H27_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H28_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H30_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H34_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
H36_sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$RAW_ROOT"castaneus_bam.txt

for bam in `cat "$RAW_ROOT"castaneus_bam.txt`
do
wget -P "$RAW_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/"$bam"
done



