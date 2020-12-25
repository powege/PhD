#!/bin/bash

# Set working directory
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/

# create file of  bam files

# Mus m musculus
# 22 individuals from 3 populations
echo '396_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
413_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
416_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
424_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
435_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
444_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17072_AL1.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17073_AL16.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17074_AL19.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17075_AL33.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17076_AL38.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17077_AL40.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17078_AL41.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570284_17079_AL42.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17122_CR12.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17123_CR13.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17124_CR16.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17125_CR23.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17126_CR25.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17127_CR29.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17128_CR46.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570285_17129_CR59.sorted.cigar.nodups.realigned.recalibrated.sorted.bam'  > "$RAW_ROOT"musculus_bam.txt

for bam in `cat "$RAW_ROOT"musculus_bam.txt`
do
wget -P "$RAW_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/"$bam"
done



