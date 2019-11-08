#!/bin/bash

### Script that downloads and unzips GWAS catalog

# set file directory
PED_ROOT=/well/lindgren/George/Data/GWAS/raw/

wget https://www.ebi.ac.uk/gwas/api/search/downloads/full \
-P "$PED_ROOT"
mv "$PED_ROOT"full "$PED_ROOT"gwas_catalog_v1.0-associations_e96_r2019-05-03.tsv

