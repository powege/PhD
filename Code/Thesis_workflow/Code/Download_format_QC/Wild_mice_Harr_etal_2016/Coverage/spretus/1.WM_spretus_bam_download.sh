#!/bin/bash

# Set working directory
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/

# create file of  bam files

# Mus spretus
# 8 individuals
echo 'SN4570286_17098_SP36.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17099_SP39.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17100_SP41.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17101_SP51.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17102_SP62.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17103_SP68.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17104_SP69.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17105_SP70.sorted.cigar.nodups.realigned.recalibrated.sorted.bam'  > "$RAW_ROOT"spretus_bam.txt

for bam in `cat "$RAW_ROOT"spretus_bam.txt`
do
wget -P "$RAW_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/"$bam"
done



