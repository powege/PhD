#!/bin/bash

### Script that downloads ClinVar data 

# set file directory
PED_ROOT=/well/lindgren/George/Data/Thesis_workflow/Data/IMPC/Raw/

# release 12
wget  ftp://ftp.ebi.ac.uk/pub/databases/impc//all-data-releases/release-12.0/results/* \
-P "$PED_ROOT"


 