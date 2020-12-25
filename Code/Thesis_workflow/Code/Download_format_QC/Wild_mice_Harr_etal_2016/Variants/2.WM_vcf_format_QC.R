### SCRIPT that formats and QCs raw vcf files from Harr et al (2016)
# QC:
  # Subset variants with alternate base == 1
  # PASS filter status

# Output:
  # gziped vcf with QCed variants by autosome and X

rm(list = ls())
graphics.off()

library(data.table)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

in.vcf <- args[1]
out.prefix <- args[2]

# in.vcf <- "/well/lindgren/George/Data/Thesis_workflow/Thesis_data/Wild_mice/Variants/Raw/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"
# out.vcf <- ""

# in.vcf <- "/well/lindgren/George/Data/Thesis_workflow/Data/Wild_mice_Harr_etal_2016/Variants/Raw/Mmd_HEL.vcf_90_recalibrated_snps_raw_indels.env.ef.SNP.MNP.vcf"
# out.vcf <- ""

### IMPORT 
all <- fread(input = in.vcf, 
             select = 1:7, 
             col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER"))

### FORAMT AND QC 

# Subset all PASS filter status
all <- all[FILTER == "PASS"] 

# Subset variants with alternate base == 1
all <- all[ALT %in% c("A", "C", "G", "T")] 

# Remove ID and QUAL
all$ID <- NA
all$QUAL <- NA

# Remove chr from CHROM
all$CHROM <- gsub("chr", "", all$CHROM)

### EXPORT

for (chr in c(1:19, "X")){
  fwrite(all[CHROM == chr], 
         file = paste0(out.prefix, chr, ".vcf.gz"), 
         compress = "gzip",
         col.names = F)
}






