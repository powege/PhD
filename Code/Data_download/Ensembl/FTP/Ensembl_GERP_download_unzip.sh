#!/bin/bash

### Script that downloads and unzips GERP constrained elements from Ensembl 
### ftp (v94, GRCh38)

# set file directory
PED_ROOT=/well/lindgren/George/Data/Ensembl/GERP/

# download and unzip mouse GERP data 
wget ftp://ftp.ensembl.org/pub/release-94/bed/ensembl-compara/70_mammals.gerp_constrained_element/gerp_constrained_elements.mus_musculus.bed.gz \
-P "$PED_ROOT"/
gunzip "$PED_ROOT"/gerp_constrained_elements.mus_musculus.bed.gz

# download and unzip human GERP data 
wget ftp://ftp.ensembl.org/pub/release-94/bed/ensembl-compara/70_mammals.gerp_constrained_element/gerp_constrained_elements.homo_sapiens.bed.gz \
-P "$PED_ROOT"/
gunzip "$PED_ROOT"/gerp_constrained_elements.homo_sapiens.bed.gz
