#!/bin/bash

### Script that downloads ClinVar data 

# set file directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/ClinVar/Raw/

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz \
-P "$PED_ROOT"

