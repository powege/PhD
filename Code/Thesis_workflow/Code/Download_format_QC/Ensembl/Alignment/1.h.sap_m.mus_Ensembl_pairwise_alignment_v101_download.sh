#!/bin/bash

########
### SET VARS
########

# set working directory for files to download into
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Raw/
cd $PED_ROOT

########
### download and unzip pairwise alignment
########

wget ftp://ftp.ensembl.org/pub/release-101/maf/ensembl-compara/pairwise_alignments/hsap_grch38.v.mmus_grcm38.lastz_net.tar.gz \
-P "$PED_ROOT"

### UNZIP 
tar xvzf "$PED_ROOT"hsap_grch38.v.mmus_grcm38.lastz_net.tar.gz

