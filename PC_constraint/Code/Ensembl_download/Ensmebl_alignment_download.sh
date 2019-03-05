#!/bin/bash

### Script that downloads and unzips ensembl allignments

# Set working directory for files to download into
cd /NGS/users/George/Workflows/PC_constraint/Paper/Data/Ensembl_allignments/

# Download Mouse Mus Caroli allignment 
#wget ftp://ftp.ensembl.org/pub/release-92/maf/ensembl-compara/pairwise_alignments/mus_musculus_GRCm38_vs_mus_caroli_CAROLI_EIJ_v1_1_lastz_net.tar.gz
wget ftp://ftp.ensembl.org/pub/release-94/maf/ensembl-compara/pairwise_alignments/mmus_grcm38.v.mcar_caroli_eij_v1.1.lastz_net.tar.gz


# Download Human Chimpanzee allignment 
#wget ftp://ftp.ensembl.org/pub/release-92/maf/ensembl-compara/pairwise_alignments/homo_sapiens_GRCh38_vs_pan_troglodytes_Pan_tro_3_0_lastz_net.tar
wget ftp://ftp.ensembl.org/pub/release-94/maf/ensembl-compara/pairwise_alignments/hsap_grch38.v.ptro_pan_tro_3.0.lastz_net.tar.gz

# unzip files
gunzip mmus_grcm38.v.mcar_caroli_eij_v1.1.lastz_net.tar.gz
gunzip hsap_grch38.v.ptro_pan_tro_3.0.lastz_net.tar.gz

# open .tar files
tar -xvf mmus_grcm38.v.mcar_caroli_eij_v1.1.lastz_net.tar
tar -xvf hsap_grch38.v.ptro_pan_tro_3.0.lastz_net.tar

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
