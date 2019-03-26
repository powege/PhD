#!/bin/bash

### Script that downloads and unzips Mouse to human alignment from Ensembl

# Set working directory for files to download into
cd /NGS/users/George/Workflows/NC_constraint/Data/Ensembl_alignments

# Download mouse to human alignment 
wget ftp://ftp.ensembl.org/pub/release-94/maf/ensembl-compara/pairwise_alignments/mmus_grcm38.v.hsap_grch38.lastz_net.tar.gz

# unzip files
gunzip mmus_grcm38.v.hsap_grch38.lastz_net.tar.gz

# open .tar files
tar -xvf mmus_grcm38.v.hsap_grch38.lastz_net.tar

## unzip files
#for FILE in mus_musculus_GRCm38_vs_mus_caroli_CAROLI_EIJ_v1_1_lastz_net/mus_musculus_GRCm38_vs_mus_caroli_CAROLI_EIJ_v1_1_lastz_net*
#do
	#echo "unzipping $FILE"
	#gunzip $FILE
	#echo "$FILE unzipped"	
#done

## unzip files
#for FILE in hsap_grch38.v.ptro_pan_tro_3.0.lastz_net/homo_sapiens_GRCh38_vs_pan_troglodytes_Pan_tro_3_0_lastz_net*
#do
	#echo "unzipping $FILE"
	#gunzip $FILE
	#echo "$FILE unzipped"	
#done
