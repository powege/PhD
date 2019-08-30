#!/bin/bash

### Script that runs Constraint_scores_canPC.R

# Run for MGP mus musculus only (ie no SPRET), and all strains. 
# Run for Human gnomAD and 1KGP, with MAF > 0.001 and 0.0001

### MOUSE
## MGP No SPRET (mus musculus only)
Rscript ~/Dropbox/GitHub_repos/PhD/Code/PC_constraint/Constraint_scores/Constraint_scores_canPC.R \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allMUSMUS_QCed_canPC_N_SNVs.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_pMu.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_nCpG.csv \
~/Dropbox/PhD/Data/MGP/coverage/Ensembl_v94_mouse_canPC_MGP_mask.csv \
~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_MGP_mus_musculus_v2.csv \
MGP

### MOUSE
## MGP All strain (including spret)
Rscript ~/Dropbox/GitHub_repos/PhD/Code/PC_constraint/Constraint_scores/Constraint_scores_canPC.R \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/MGP_v5_allSTRAINS_QCed_canPC_N_SNVs.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_pMu.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_mouse_canPC_QCpass_nCpG.csv \
~/Dropbox/PhD/Data/MGP/coverage/Ensembl_v94_mouse_canPC_MGP_mask.csv \
~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_MGP_all_strains_v2.csv \
MGP

### HUMAN
### 1KGP MAF > 0.001
Rscript ~/Dropbox/GitHub_repos/PhD/Code/PC_constraint/Constraint_scores/Constraint_scores_canPC.R \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF001.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv \
~/Dropbox/PhD/Data/1KGP/Masks/Formatted/Ensembl_v94_human_canPC_1KGP_mask.csv \
~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF001_v2.csv \
1KGP

### HUMAN
### 1KGP MAF > 0.0005
Rscript ~/Dropbox/GitHub_repos/PhD/Code/PC_constraint/Constraint_scores/Constraint_scores_canPC.R \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF0005.csv  \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv \
~/Dropbox/PhD/Data/1KGP/Masks/Formatted/Ensembl_v94_human_canPC_1KGP_mask.csv \
~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF0005_v2.csv \
1KGP

### HUMAN
### 1KGP MAF > 0.0001
Rscript ~/Dropbox/GitHub_repos/PhD/Code/PC_constraint/Constraint_scores/Constraint_scores_canPC.R \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/1KGP_phase3_QCed_canPC_N_SNVs_MAF0001.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv \
~/Dropbox/PhD/Data/1KGP/Masks/Formatted/Ensembl_v94_human_canPC_1KGP_mask.csv \
~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF0001_v2.csv \
1KGP

### HUMAN
### gnomAD MAF > 0.001
Rscript ~/Dropbox/GitHub_repos/PhD/Code/PC_constraint/Constraint_scores/Constraint_scores_canPC.R \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/gnomAD_v2.1.1_QCed_canPC_N_SNVs_MAF001.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_pMu.csv \
~/Dropbox/PhD/Data/PC_constraint/Model_variables/Ensembl_v94_human_canPC_QCpass_nCpG.csv \
~/Dropbox/PhD/Data/gnomAD/coverage/Ensembl_v94_human_canPC_gnomad_mask.csv \
~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_gnomAD_MAF001_v2.csv \
gnomAD


