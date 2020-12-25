#!/bin/bash

### LOAD MODULES
module load Kent_tools/401-gompi-2019b

### SET VARS
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Methylation/Raw/ # set working directory for files to download into

### download methylation big bed files
wget http://ftp.ensembl.org/pub/data_files/mus_musculus/GRCm38/dna_methylation_feature/ES_5mC_Publication_Stadler2011_PMID22170606/ES_5mC_Stadler2011_PMID22170606.bb \
-P "$PED_ROOT"
wget http://ftp.ensembl.org/pub/data_files/mus_musculus/GRCm38/dna_methylation_feature/NPC_5mC_Publication_Stadler2011_PMID22170606/NPC_5mC_Stadler2011_PMID22170606.bb \
-P "$PED_ROOT"

### convert bigBed to bed
bigBedToBed "$PED_ROOT"ES_5mC_Stadler2011_PMID22170606.bb "$PED_ROOT"ES_5mC_Stadler2011_PMID22170606.bed
bigBedToBed "$PED_ROOT"NPC_5mC_Stadler2011_PMID22170606.bb "$PED_ROOT"NPC_5mC_Stadler2011_PMID22170606.bed

### zip
gzip "$PED_ROOT"ES_5mC_Stadler2011_PMID22170606.bed
gzip "$PED_ROOT"NPC_5mC_Stadler2011_PMID22170606.bed



