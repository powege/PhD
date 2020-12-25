#!/bin/bash

# create files with classical strains 
# download bam files for all classical strains

# Set working directory
RAW_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/MGP/Coverage/Raw/

# create file of classical strains in vcf
echo '129P2_OlaHsd
129S1_SvImJ
129S5SvEvBrd
AKR_J
A_J
BALB_cJ
BTBR_T__Itpr3tf_J
BUB_BnJ
C3H_HeH
C3H_HeJ
C57BL_10J
C57BL_6NJ
C57BR_cdJ
C57L_J
C58_J
CBA_J
DBA_1J
DBA_2J
FVB_NJ
I_LnJ
KK_HiJ
LP_J
NOD_ShiLtJ
NZB_B1NJ
NZO_HlLtJ
NZW_LacJ
RF_J
SEA_GnJ
ST_bJ' > "$RAW_ROOT"classicalSTRAINS.txt

# download .bam files from ftp

# classical
for bam in `cat "$RAW_ROOT"classicalSTRAINS.txt`
do
wget -P "$RAW_ROOT" ftp://ftp-mouse.sanger.ac.uk/current_bams/"$bam".bam
done


