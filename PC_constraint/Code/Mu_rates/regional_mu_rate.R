# Script to estimate regional mutation rate by counting the number of SNVs 
# within 1Mbp upstream and downstream of the transcript start and end coordinates, 
# and dividing by 2,000,000.

rm(list=ls())
graphics.off()

library(data.table)


# Function that calculates regional mu rate by dividing all SNVs in region by the total bases. 
all.region <- function(ant, reg, chr){
  
  # subset by chromosome
  reg <- subset(reg, reg$chromosome_name == chr)
  
  # add 1mb to start and end coordinates
  reg$start_position_mb <- reg$start_position-1000000
  reg$end_position_mb <- reg$end_position+1000000
  
    # remove functional variants
    # ant <- subset(ant, !(ant$Consequence %like% "missense_variant") | 
    #                 !(ant$Consequence %like% "stop_gained") |
    #                 !(ant$Consequence %like% "start_lost") |
    #                 !(ant$Consequence %like% "stop_lost") |
    #                 !(ant$Consequence %like% "splice_donor_variant") |
    #                 !(ant$Consequence %like% "splice_acceptor_variant"))
  
    reg$n_SNP <- NA
    for (i in 1:nrow(reg)){
      reg$n_SNP[i] <- nrow(ant[which(ant$V2 > reg$start_position_mb[i] & ant$V2 < reg$end_position_mb[i]),])
      print(i)
    }
    
    reg$n_base <- reg$end_position_mb - reg$start_position_mb
    reg$r_mu <- reg$n_SNP/reg$n_base

  return(reg)
}

#### RUN for MOUSE

M_reg <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/M_gene_POS.csv")
CHR <- c(1:19, "X")
M_out <- data.frame()
for (i in 1:length(CHR)){
  print(paste("PROCESSING CHROMOSOME", CHR[i], sep = " "))
  M_ant.file <- paste0("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/M_MGP_QCed_VEP_all_chr", CHR[i], ".txt")
  M_ant <- fread(M_ant.file)
  tmp <- all.region(ant = M_ant,
                    reg = M_reg,
                    chr = CHR[i])
  M_out <- rbind(M_out, tmp)
  print(paste("CHORMOSOME", CHR[i], "DONE!", sep = " "))
}

fwrite(M_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_regional_mu_rates.csv")


#### RUN for HUMAN

H_reg <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_gene_POS.csv")
CHR <- c(1:22, "X")
H_out <- data.frame()
for (i in 1:length(CHR)){
  print(paste("PROCESSING CHROMOSOME", CHR[i], sep = " "))
  H_ant.file <- paste0("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/H_1000GP_QCed_VEP_all_chr", CHR[i], ".txt")
  H_ant <- fread(H_ant.file, fill = T)
  tmp <- all.region(ant = H_ant,
                    reg = H_reg,
                    chr = CHR[i])
  H_out <- rbind(H_out, tmp)
  print(paste("CHORMOSOME", CHR[i], "DONE!", sep = " "))
}
fwrite(H_out, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_regional_mu_rates.csv")


