#!/bin/bash

### Script that downloads and unzips human-mouse alignment from Ensembl 
### ftp (v94, GRCh38)

# set file directory
PED_ROOT=/well/lindgren/George/Data/Ensembl/Alignment/Raw/

# download and unzip data 
wget ftp://ftp.ensembl.org/pub/release-94/maf/ensembl-compara/pairwise_alignments/mmus_grcm38.v.hsap_grch38.lastz_net.tar.gz \
-P "$PED_ROOT"/
tar xvzf "$PED_ROOT"/mmus_grcm38.v.hsap_grch38.lastz_net.tar.gz