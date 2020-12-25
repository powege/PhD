rm(list=ls())
graphics.off()

library(data.table)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

# set args variables
in.CG.file <- args[1]
in.ES.file <- args[2]
in.NPC.file <- args[3]
out.file <- args[4]

# in.CG.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/m.mus_reference_Ensembl_GRC38_CG.bed.gz"
# in.ES.file <- "~/Dropbox/PhD/Data/Ensembl/Methylation/raw/ES_5mC_Stadler2011_PMID22170606.bed.gz"
# in.NPC.file <- "~/Dropbox/PhD/Data/Ensembl/Methylation/raw/NPC_5mC_Stadler2011_PMID22170606.bed.gz"
# out.file <- "~/Dropbox/PhD/Data/Ensembl/Methylation/formatted/m.mus_GRC38_Ensembl_CG_meth_NPC_ES.bed.gz"

### IMPORT

CG <- fread(paste0("gunzip -cq ", in.CG.file))
ES <- fread(paste0("gunzip -cq ", in.ES.file))
NPC <- fread(paste0("gunzip -cq ", in.NPC.file))

### FORMAT 

colnames(CG) <- c("chromosome", "start", "end")
colnames(ES) <- c("chromosome", "start", "end", "context_coverage", "ES_percent_methylated", "strand")
colnames(NPC) <- c("chromosome", "start", "end", "context_coverage", "NPC_percent_methylated", "strand")

ES$ES_percent_methylated <- ES$ES_percent_methylated / 10
NPC$NPC_percent_methylated <- NPC$NPC_percent_methylated / 10

# split context_coverage column
ES[, c("ES_context", "ES_coverage") := tstrsplit(context_coverage, "/", fixed=TRUE)][, context_coverage := NULL]
NPC[, c("NPC_context", "NPC_coverage") := tstrsplit(context_coverage, "/", fixed=TRUE)][, context_coverage := NULL]

ES$chromosome <- gsub("chr", "", ES$chromosome)
NPC$chromosome <- gsub("chr", "", NPC$chromosome)

ES <- ES[chromosome %in% as.character(1:19)]
NPC <- NPC[chromosome %in% as.character(1:19)]
CG$chromosome <- as.character(CG$chromosome)

CG <- CG[order(chromosome, start),]
ES <- ES[order(chromosome, start),]
NPC <- NPC[order(chromosome, start),]

# correct for mismatch in alignment
ES$start <- ES$end
ES$end <- ES$start + 1
NPC$start <- NPC$end
NPC$end <- NPC$start + 1

# merge
output <- ES[CG, on = c("chromosome", "start", "end")]
output <- NPC[output, on = c("chromosome", "start", "end", "strand")]
output <- output[!duplicated(output[,1:3]),]

# average coverage and % methylated between ES and NPC
output$coverage <- as.numeric(output$NPC_coverage) + as.numeric(output$ES_coverage)
output$percent_methylated <-  (( (as.numeric(output$ES_coverage) * (output$ES_percent_methylated / 100)) + 
                                 (as.numeric(output$NPC_coverage) * (output$NPC_percent_methylated / 100)) ) / output$coverage ) *100
                               
output <- output[,c("chromosome",
                    "start",
                    "end",
                    "NPC_percent_methylated",
                    "NPC_coverage",
                    "ES_percent_methylated",
                    "ES_coverage",
                    "coverage",
                    "percent_methylated")]

output[is.na(output)] <- 0 # replace NA with 0
gw_coverage <- (nrow(output[coverage != 0]) / nrow(output)) *100 # % CGs covered by ES and NPC
ES_NPC_cor <- cor.test(output$NPC_percent_methylated, output$ES_percent_methylated)$estimate # correlation between ES and NPC % methylated

### EXPORT

fwrite(output, out.file, compress = "gzip")


#####

# str(output)
# summary(output)
# 
# out_all <- output
# out_all[is.na(out_all)] <- 0
# 
# hist(output$percent_methylated)
# hist(out_all$percent_methylated)
# 
# cor.test(output$NPC_percent_methylated, output$ES_percent_methylated)
# cor.test(out_all$NPC_percent_methylated, out_all$ES_percent_methylated)


# The distribution of mean methylation values for CG dinucleotides across the mouse genome, 
# sampled from embryonic stem cells.  
# We divided the genome into 3 levels (low methylation, missing or < 0.2; medium, 0.2-0.6; 
# and high, > 0.6) and computed all ensuing metrics based on these categories.

