### SCRIPT that formats and QCs raw vcf files from MGP
# QC:
  # Subset variants with alternate base == 1
  # PASS filter status
  # homozygous in one or more strain

# Output:
  # gziped vcf with QCed variants for all strains
  # gziped vcf with QCed variants for all clasical strains
  # gziped vcf with QCed variants for all wild-derived strains


rm(list = ls())
graphics.off()

library(data.table)


### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

in.vcf <- args[1]
out.vcf.all <- args[2]
out.vcf.classical <- args[3]
out.vcf.wild <- args[4]

# in.vcf <- "/well/lindgren/George/Data/Thesis_workflow/Thesis_data/MGP/Variants/Raw/mgp.v5.merged.snps_all.dbSNP142_chr17.vcf.gz"
# out.vcf.all <- ""
# out.vcf.classical <- ""
# out.vcf.wild <- ""


### SET VARS

vcf_cols_all <- c("CHROM",
                   "POS",
                   "ID",
                   "REF",
                   "ALT",
                   "QUAL",
                   "FILTER",
                   "INFO",
                   "FORMAT",
                   "129P2_OlaHsd",
                   "129S1_SvImJ",
                   "129S5SvEvBrd",
                   "AKR_J",
                   "A_J",
                   "BALB_cJ",
                   "BTBR_T+_Itpr3tf_J",
                   "BUB_BnJ",
                   "C3H_HeH",
                   "C3H_HeJ",
                   "C57BL_10J",
                   "C57BL_6NJ",
                   "C57BR_cdJ",
                   "C57L_J",
                   "C58_J",
                   "CAST_EiJ",
                   "CBA_J",
                   "DBA_1J",
                   "DBA_2J",
                   "FVB_NJ",
                   "I_LnJ",
                   "KK_HiJ",
                   "LEWES_EiJ",
                   "LP_J",
                   "MOLF_EiJ",
                   "NOD_ShiLtJ",
                   "NZB_B1NJ",
                   "NZO_HlLtJ",
                   "NZW_LacJ",
                   "PWK_PhJ",
                   "RF_J",
                   "SEA_GnJ",
                   "SPRET_EiJ",
                   "ST_bJ",
                   "WSB_EiJ",
                   "ZALENDE_EiJ")
vcf_cols <- c("CHROM",
              "POS",
              "ID",
              "REF",
              "ALT",
              "QUAL",
              "FILTER",
              "INFO",
              "FORMAT")
wild_strains <- c("CAST_EiJ",
                  "WSB_EiJ",
                  "PWK_PhJ",
                  "MOLF_EiJ",
                  "SPRET_EiJ",
                  "LEWES_EiJ",
                  "ZALENDE_EiJ")
classical_strains <- vcf_cols_all[!vcf_cols_all %in% c(vcf_cols, wild_strains)]
wild_strain_cols <- c(vcf_cols, wild_strains)
classical_strain_cols <- vcf_cols_all[!vcf_cols_all %in% wild_strains]


### IMPORT 

all <- fread(paste0("gunzip -cq ", in.vcf), fill = T, sep = "\t")


### FORAMT AND QC 

colnames(all) <- vcf_cols_all

# Subset all PASS filter status
all <- all[FILTER == "PASS"] 

# Subset variants with alternate base == 1
all <- all[ALT %in% c("A", "C", "G", "T")] 

# subset vcf columns with genotype 1/1 in one or more wild strain
wild <- all[apply(all[, ..wild_strains], 1, function(i) any(grepl("1/1", i))), 1:7] 

# subset vcf columns with genotype 1/1 in one or more classical strain
classical <- all[apply(all[, ..classical_strains], 1, function(i) any(grepl("1/1", i))), 1:7] 

# subset vcf columns for all
all <- all[,1:7]

### EXPORT

fwrite(all, out.vcf.all, compress = "gzip")
fwrite(wild, out.vcf.wild, compress = "gzip")
fwrite(classical, out.vcf.classical, compress = "gzip")





