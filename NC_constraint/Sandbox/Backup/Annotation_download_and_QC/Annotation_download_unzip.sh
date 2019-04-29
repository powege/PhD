#!/bin/bash

### Script that downloads and unzips human and mouse genome-wide annotations from
### Ensembl ftp (v94, GRCh38), including: genes, Regulatory Build, and motifs. 

# set file directory
PED_ROOT=/well/lindgren/George/Data/Ensembl/

# download and unzip human gene annotation
wget ftp://ftp.ensembl.org/pub/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh38.94.gff3.gz \
-P "$PED_ROOT"/Annotation/
gunzip "$PED_ROOT"/Annotation/Homo_sapiens.GRCh38.94.gff3.gz

# download and unzip mouse gene annotation
wget ftp://ftp.ensembl.org/pub/release-94/gff3/mus_musculus/Mus_musculus.GRCm38.94.gff3.gz \
-P "$PED_ROOT"/Annotation/
gunzip "$PED_ROOT"/Annotation/Mus_musculus.GRCm38.94.gff3.gz

# download and unzip human Regulatory Build annotation
wget ftp://ftp.ensembl.org/pub/release-94/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz \
-P "$PED_ROOT"/Annotation/
gunzip "$PED_ROOT"/Annotation/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz

# download and unzip mouse Regulatory Build annotation
wget ftp://ftp.ensembl.org/pub/release-94/regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz \
-P "$PED_ROOT"/Annotation/
gunzip "$PED_ROOT"/Annotation/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz

# download and unzip human motif annotation
wget ftp://ftp.ensembl.org/pub/current_regulation/homo_sapiens/MotifFeatures/Homo_sapiens.GRCh38.motif_features.gff.gz \
-P "$PED_ROOT"/Annotation/
gunzip "$PED_ROOT"/Annotation/Homo_sapiens.GRCh38.motif_features.gff.gz

# download and unzip mouse motif annotation
wget ftp://ftp.ensembl.org/pub/current_regulation/mus_musculus/MotifFeatures/Mus_musculus.GRCm38.motif_features.gff.gz \
-P "$PED_ROOT"/Annotation/
gunzip "$PED_ROOT"/Annotation/Mus_musculus.GRCm38.motif_features.gff.gz


