### Table 1 – Comparison of total genes and SNV counts in the Mouse Genomes Project dataset (mouse), and the 1000 Genomes Project dataset (human) at different minor allele frequency (MAF) thresholds. 

rm(list = ls())
graphics.off()

library(gridExtra)
library(data.table)
library(plyr)

### IMPORT
H <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/1000GP_canPC_SNV_totals.txt")
M <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/SNV_counts/MGP_canPC_SNV_totals.txt")


### FORMAT

H <- H[H$MAF_threshold %in% c(0.0001, 0.0005, 0.001),]
H$MAF_threshold <- c("0.0001", "0.0005", "0.001")
M$MAF_threshold <- "NA"
H$Species <- "Human"
M$Species <- "Mouse"

colnames(H) <- c("MAF > threshold","N genes", "N synonymous", "N missense", "N nonsense",  
                 "N splice acceptor/donor", "N intron", "Species")
colnames(M) <- c("N genes", "N synonymous", "N missense", "N nonsense",  
                 "N splice acceptor/donor", "N intron", "MAF > threshold", "Species")

out <- rbind.fill(M, H)
out <- out[,c("Species", "MAF > threshold", "N genes", "N synonymous", "N missense", "N nonsense",  
              "N splice acceptor/donor", "N intron")]

png("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Tables/Table_1.pdf", height = 30*nrow(out), width = 90*ncol(out))
grid.table(out)
dev.off()
 