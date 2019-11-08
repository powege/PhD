#!/bin/bash

### Script that downloads and unzips ClinVar data 

# set file directory
PED_ROOT=/well/lindgren/George/Data/ClinVar/raw/

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz \
-P "$PED_ROOT"
gunzip "$PED_ROOT"variant_summary.txt.gz



