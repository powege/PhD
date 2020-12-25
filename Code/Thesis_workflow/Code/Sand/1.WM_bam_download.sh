#!/bin/bash

# Script that downloads .bam files for all mus musculus domesticus individuals

# Set working directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Coverage/Raw/

# create file of wild-derived strain bam files in vcf

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
H36_sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$PED_ROOT"castaneus.txt

# Mus m domesticus
# 27 individuals (24 individuals from 3 populations of Mus m domesticus, and 3 individuals of Mus m helgolandicus)
# Mus m helgolandicus is being treated as Mus m domesticus for analyses due to molecular proximity, as done by Harr et al (2016).
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
TP81B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN7640133_6820_HG06.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN7640133_6821_HG08.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN7640133_6822_HG13.sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$PED_ROOT"domesticus.txt

### NB: Mus mus helgolandicus is being treated as Mus mus domesticus group for analyses due to molecular proximity, 
# as done by Harr et al (2016).
# echo 'SN7640133_6820_HG06.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
# SN7640133_6821_HG08.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
# SN7640133_6822_HG13.sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$PED_ROOT"helgolandicus.txt

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
SN4570285_17129_CR59.sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$PED_ROOT"musculus.txt

# Mus spretus
# 8 individuals
echo 'SN4570286_17098_SP36.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17099_SP39.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17100_SP41.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17101_SP51.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17102_SP62.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17103_SP68.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17104_SP69.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
SN4570286_17105_SP70.sorted.cigar.nodups.realigned.recalibrated.sorted.bam' > "$PED_ROOT"spretus.txt


# download files

# M m castaneus

wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H12_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H14_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H15_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H24_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H26_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H27_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H28_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H30_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H34_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/H36_sorted.cigar.nodups.realigned.recalibrated.sorted.bam

# M m domesticus

wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/14_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/15B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/16B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/18B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/AH15_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/AH23_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/B2C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/C1_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/E1_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/F1B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/JR11_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/JR15_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/JR2-F1C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/JR5-F1C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/JR7-F1C_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/JR8-F1A_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP121B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP17-2_ChrAdded_ordered.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP1_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP3-92_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP4a_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP51D_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP7-10F1A2_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/TP81B_sorted.cigar.nodups.realigned.recalibrated.sorted.bam

# M m helgolandicus

wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_helgolandicus/genomes_bam/SN7640133_6820_HG06.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_helgolandicus/genomes_bam/SN7640133_6821_HG08.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_helgolandicus/genomes_bam/SN7640133_6822_HG13.sorted.cigar.nodups.realigned.recalibrated.sorted.bam

# M m musculus

wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/396_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/413_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/416_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/424_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/435_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/444_sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17072_AL1.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17073_AL16.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17074_AL19.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17075_AL33.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17076_AL38.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17077_AL40.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17078_AL41.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570284_17079_AL42.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17122_CR12.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17123_CR13.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17124_CR16.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17125_CR23.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17126_CR25.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17127_CR29.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17128_CR46.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/SN4570285_17129_CR59.sorted.cigar.nodups.realigned.recalibrated.sorted.bam

# M spretus

wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17098_SP36.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17099_SP39.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17100_SP41.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17101_SP51.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17102_SP62.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17103_SP68.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17104_SP69.sorted.cigar.nodups.realigned.recalibrated.sorted.bam
wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/SN4570286_17105_SP70.sorted.cigar.nodups.realigned.recalibrated.sorted.bam


##### 

# curl does not work for http!!! 
#
# # M m castaneus
# curl -l http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/ > "$PED_ROOT"castaneus.txt
# grep ".bam$" "$PED_ROOT"castaneus.txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"castaneus.txt
# for bam in `cat "$PED_ROOT"castaneus.txt`
# do
# wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_castaneus/genomes_bam/"$bam"
# done
# 
# # M m domesticus
# curl -l http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/ > "$PED_ROOT"domesticus.txt
# grep ".bam$" "$PED_ROOT"domesticus.txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"domesticus.txt
# for bam in `cat "$PED_ROOT"domesticus.txt`
# do
# wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_domesticus/genomes_bam/"$bam"
# done
# 
# # M m helgolandicus
# curl -l http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_helgolandicus/genomes_bam/ > "$PED_ROOT"helgolandicus.txt
# grep ".bam$" "$PED_ROOT"helgolandicus.txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"helgolandicus.txt
# for bam in `cat "$PED_ROOT"helgolandicus.txt`
# do
# wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_helgolandicus/genomes_bam/"$bam"
# done
# 
# # M m musculus
# curl -l http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/ > "$PED_ROOT"musculus.txt
# grep ".bam$" "$PED_ROOT"musculus.txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"musculus.txt
# for bam in `cat "$PED_ROOT"musculus.txt`
# do
# wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_m_musculus/genomes_bam/"$bam"
# done
# 
# # M m spretus
# curl -l http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/ > "$PED_ROOT"spretus.txt
# grep ".bam$" "$PED_ROOT"spretus.txt > "$PED_ROOT"tmp && mv "$PED_ROOT"tmp "$PED_ROOT"spretus.txt
# for bam in `cat "$PED_ROOT"spretus.txt`
# do
# wget -P "$PED_ROOT" http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/m_spretus/genomes_bam/"$bam"
# done


# # list files
# 
# cd "$PED_ROOT"
# ls . > "$PED_ROOT"FILES.txt
# sed -i '/FILES.txt/d' "$PED_ROOT"FILES.txt