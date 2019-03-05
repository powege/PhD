# Script that estimates mutant and ancestral sequences

rm(list =ls())
graphics.off()

library(data.table)
library(dplyr)


### FUNCTION that estimates intergenic mutant and ancestral sequences and SNVs
ancestral_mutant <- function(allignment_in_path,
                             annotation_in_path,
                             allignment_in_file,
                             annotation_in_file,
                             seq_out_path,
                             SNV_out_path,
                             seq_out_file,
                             SNV_out_file,
                             exon.pos, 
                             chr,
                             colnames_ant){
  
  ### IMPORT
  
  # sequence allignment
  allignment_file_path <- paste0(allignment_in_path, allignment_in_file, chr, ".txt")
  allignment <- fread(allignment_file_path)
  # annotations 
  ant_file_path <- paste0(annotation_in_path, annotation_in_file, chr, ".txt")
  ant <- fread(ant_file_path, fill = T)
  
  
  ### FORMAT
  
  # colnames
  colnames(allignment) <- c("POS", "REF", "SPEC")
  colnames(ant) <- colnames_ant
  
  # # convert NA to "-"
  # allignment$REF[is.na(allignment$REF)] <- "-"
  # allignment$SPEC[is.na(allignment$SPEC)] <- "-"
  
  # convert "" to "-"
  allignment$REF[allignment$REF == ""] <- "-"
  allignment$SPEC[allignment$SPEC == ""] <- "-"
  
  ## cut exons from df sequence 
  exon.pos <- subset(exon.pos, exon.pos$chromosome_name == chr)
  
  # creat list of exon sequences (vectors)
  exon.pos.list <- list()
  for (i in 1:nrow(exon.pos)){
    exon.pos.list[[i]] <- seq(from = exon.pos$exon_chrom_start[i], to = exon.pos$exon_chrom_end[i])
  }
  exon.pos.list <- unlist(exon.pos.list)
  
  # cut all positions within exons
  seq.exon.cut <- allignment$POS[!allignment$POS %in% exon.pos.list]
  seq.exon.cut <- data.table(POS = seq.exon.cut)
  allignment <- allignment[seq.exon.cut, on = "POS", nomatch = 0]
  rm(seq.exon.cut)
  
  ## cut exon POS from annotations
  ant <- ant[, c("POS", "REF", "ALT"), with = FALSE]
  # test1 <- allignment[allignment$POS %in% ant$POS,]
  # test2 <- ant[ant$POS %in% test1$POS,]
  # test3 <- cbind(test1, test2)
  # test <- ant[allignment, on = c("POS", "REF"), nomatch=0]
  ant <- ant[allignment, on = c("POS", "REF"), nomatch=0]
  
  
  ### ANCESTRAL sequence
  
  # Assume “ancestral” base is REF, or ALT if ALT is shared with related species. 
  ant$ANCESTRAL <- ifelse(ant$ALT == ant$SPEC, ant$ALT, ant$REF) # vectorisation
  ind <- which(allignment$POS %in% ant$POS)
  allignment$ANCESTRAL <- replace(allignment$REF, ind, ant$ANCESTRAL)
  
  ### MUTANT sequence
  
  # Assume “mutant” base is ALT, or REF if ALT is shared with related species.
  ant$MUTANT <- ifelse(ant$ALT == ant$SPEC, ant$REF, ant$ALT) # vectorisation
  # ind <- which(allignment$POS %in% ant$POS)
  allignment$MUTANT <- replace(allignment$REF, ind, ant$MUTANT)
  
  # Subset output columns
  SEQ <- allignment[,c("POS", "ANCESTRAL", "MUTANT")]
  SNV <- ant[,c("POS", "ANCESTRAL", "MUTANT")]
  
  # remove all rows where ANCESTRAL == "-" | MUTANT == "-"
  rm.id <- which((SEQ$ANCESTRAL == "-" | SEQ$MUTANT == "-"))
  if (length(rm.id) != 0){
    SEQ <- SEQ[-rm.id,]
  }
  rm.id <- which((SNV$ANCESTRAL == "-" | SNV$MUTANT == "-"))
  if (length(rm.id) != 0){
    SNV <- SNV[-rm.id,]
  }
  
  # split into character strings of sequential POS
  # ANCESTRAL <- split(allignment$ANCESTRAL, cumsum(c(TRUE, diff(allignment$POS)!=1)))
  # MUTANT <- split(allignment$MUTANT, cumsum(c(TRUE, diff(allignment$POS)!=1)))
  # ANCESTRAL <- lapply(X=ANCESTRAL, FUN=paste, collapse = '')
  
  ### OUTPUT
  SEQ_file_path <- paste0(seq_out_path, seq_out_file, chr, ".txt")
  SNV_file_path <- paste0(SNV_out_path, SNV_out_file, chr, ".txt")
  fwrite(SEQ, SEQ_file_path, col.names = F, sep = '\t')
  fwrite(SNV, SNV_file_path, col.names = F, sep = '\t')
  
  print(paste("chromosome", chr, "dome!", sep = " "))
  
}


### RUN FOR MOUSE
# M_exon.pos <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/M_exon_POS.csv")
# for (CHR in c(1:19, "X")){
# ancestral_mutant(allignment_in_path = "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_alignments/",
#                  annotation_in_path = "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/",
#                  allignment_in_file = "Mouse_Muscaroli_alignment_chr",
#                  annotation_in_file = "M_MGP_QCed_VEP_all_chr",
#                  seq_out_path = "/Volumes/Untitled/PC_constraint/Paper/Data/Ancestral_mutant/",
#                  SNV_out_path = "/Volumes/Untitled/PC_constraint/Paper/Data/Ancestral_mutant/",
#                  seq_out_file = "M_ancestral_mutant_SEQ_chr",
#                  SNV_out_file = "M_ancestral_mutant_SNP_chr",
#                  exon.pos = M_exon.pos,
#                  chr = CHR,
#                  colnames_ant =  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "Gene", "Feature", "Feature_type", 
#                                    "Consequence", "IMPACT", "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS"))
# }

### RUN FOR HUMAN
H_exon.pos <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Ensembl/H_exon_POS.csv")
for (CHR in c(1:22, "X")){
  ancestral_mutant(allignment_in_path = "/Volumes/HarEx/PC_constraint/Paper/Data/Ensembl_alignments/",
                   annotation_in_path = "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/",
                   allignment_in_file = "Human_chimpanzee_alignment_chr",
                   annotation_in_file = "H_1000GP_QCed_VEP_all_chr",
                   seq_out_path = "/Volumes/HarEx/PC_constraint/Paper/Data/Ancestral_mutant/",
                   SNV_out_path = "/Volumes/HarEx/PC_constraint/Paper/Data/Ancestral_mutant/",
                   seq_out_file = "H_ancestral_mutant_SEQ_chr",
                   SNV_out_file = "H_ancestral_mutant_SNP_chr",
                   exon.pos = H_exon.pos,
                   chr = CHR,
                   colnames_ant =  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC", "Gene", "Feature", "Feature_type", 
                                     "Consequence", "IMPACT", "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS", "AF", 
                                     "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF"))
}

#####

### test 
# allignment_in_path = "/Volumes/HarEx/PC_constraint/Paper/Data/Ensembl_alignments/"
# annotation_in_path = "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/"
# allignment_in_file = "Human_chimpanzee_alignment_chr"
# annotation_in_file = "H_1000GP_QCed_VEP_all_chr"
# seq_out_path = "/Volumes/HarEx/PC_constraint/Paper/Data/Ancestral_mutant/"
# SNV_out_path = "/Volumes/HarEx/PC_constraint/Paper/Data/Ancestral_mutant/"
# seq_out_file = "H_ancestral_mutant_SEQ_chr"
# SNV_out_file = "H_ancestral_mutant_SNP_chr"
# exon.pos = H_exon.pos
# chr = 22
# colnames_ant =  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "AC", "Gene", "Feature", "Feature_type",
#                   "Consequence", "IMPACT", "SYMBOL", "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "CCDS", "AF", 
#                   "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF")

