#!/bin/bash

# Set working directory
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/

# create file of  bam files

# Mus m domesticus
# 24 individuals from 3 populations of Mus m domesticus
echo '14_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
15B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
16B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
18B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
AH15_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
AH23_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
B2C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
C1_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
E1_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
F1B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
JR11_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
JR15_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
JR2-F1C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
JR5-F1C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
JR7-F1C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
JR8-F1A_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP121B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP17-2_ChrAdded_ordered.bam
TP1_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP3-92_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP4a_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP51D_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP7-10F1A2_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
TP81B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$RAW_ROOT"domesticus_bam.txt

# Mus m helgolandicus
# 3 individuals 
echo 'SN7640133_6820_HG06.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN7640133_6821_HG08.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN7640133_6822_HG13.sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$RAW_ROOT"helgolandicus_bam.txt

# Mus m domesticus helgolandicus
# 27 individuals (24 individuals from 3 populations of Mus m domesticus, and 3 individuals of Mus m helgolandicus)
# Mus m helgolandicus is being treated as Mus m domesticus for analyses due to molecular proximity, as done by Harr et al (2016).
cat "$RAW_ROOT"domesticus_bam.txt "$RAW_ROOT"helgolandicus_bam.txt > "$RAW_ROOT"domesticus_helgolandicus_bam.txt

# download domesticus
for bam in `cat "$RAW_ROOT"domesticus_bam.txt`
do
wget -P "$RAW_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/"$bam"
done

# download helgolandicus
for bam in `cat "$RAW_ROOT"helgolandicus_bam.txt`
do
wget -P "$RAW_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_helgolandicus/genomes_bam/"$bam"
done


