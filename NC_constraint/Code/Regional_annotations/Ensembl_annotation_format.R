rm(list=ls())
graphics.off()

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)

### FUNCTIONS

# sub <- subset(exon_pos, exon_pos$ensembl_transcript_id == "ENSMUST00000027736")
### FUNCTION that infers intron POS from exon POS
get_intron_POS <- function(sub){
  if (nrow(sub) == 1){
    out <- data.frame(chromosome = sub$chromosome_name,
                      start = NA,
                      end = NA)
  }
  if (nrow(sub) > 1){
    sub <- sub[order(sub$exon_chrom_start),] 
    start <- sub$exon_chrom_end[1:(length(sub$exon_chrom_end))-1]
    end <- sub$exon_chrom_start[2:(length(sub$exon_chrom_end))]
    chrom <- rep(sub$chromosome_name[1], length(start))
    out <- data.frame(chromosome = chrom,
                      start = start,
                      end = end)
  }
  # dif <- sub$exon_chrom_end - sub$exon_chrom_start
  # len <- sum(dif) + length(dif)
  # out <- c(sub$ensembl_transcript_id[1], len)
  # names(out) <- c("ensembl_transcript_id", "intron_length")
  return(out)
}


### IMPORT Ensembl data
reg_build <-  fread("~/Dropbox/PhD/Data/Ensembl/Regulatory_build/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff")
exon_pos <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/M_exon_POS.csv")
UTR_pos <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/M_UTR_POS.csv")

### FORMAT .gff 
colnames(reg_build) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
reg_build <- 
  reg_build %>% separate(attribute, c("ID", "bound_end","bound_start", "description", "feature_type"), sep = ";")
reg_build$ID <- gsub("ID=", "", reg_build$ID)
reg_build$bound_end <- gsub("bound_end=", "", reg_build$bound_end)
reg_build$bound_start <- gsub("bound_start=", "", reg_build$bound_start)
reg_build$description <- gsub("description=", "", reg_build$description)
reg_build$feature_type <- gsub("feature_type=", "", reg_build$feature_type)

# table feature types
# table(reg_build$feature_type)

# subset promotor POS
promoter_pos <- subset(reg_build, reg_build$feature_type == "Promoter")
promoter_pos <- promoter_pos[,c("chromosome", "start", "end")]
promoter_pos <- promoter_pos[complete.cases(promoter_pos),]

# subset promotor flanking POS
promoter_flank_pos <- subset(reg_build, reg_build$feature_type == "Promoter Flanking Region")
promoter_flank_pos <- promoter_flank_pos[,c("chromosome", "start", "end")]
promoter_flank_pos <- promoter_flank_pos[complete.cases(promoter_flank_pos),]

# subset enhancer POS
enhancer_pos <- subset(reg_build, reg_build$feature_type == "Enhancer")
enhancer_pos <- enhancer_pos[,c("chromosome", "start", "end")]
enhancer_pos <- enhancer_pos[complete.cases(enhancer_pos),]

# subset TF binding site POS
TF_binding_pos <- subset(reg_build, reg_build$feature_type == "TF binding site" | reg_build$feature_type == "CTCF Binding Site")
TF_binding_pos <- TF_binding_pos[,c("chromosome", "start", "end")]
TF_binding_pos <- TF_binding_pos[complete.cases(TF_binding_pos),]

# subset open chromatin POS
open_chromatin_pos <- subset(reg_build, reg_build$feature_type == "Open chromatin")
open_chromatin_pos <- open_chromatin_pos[,c("chromosome", "start", "end")]
open_chromatin_pos <- open_chromatin_pos[complete.cases(open_chromatin_pos),]

# subset PC exon POS
PC_exon_pos <- subset(exon_pos, exon_pos$gene_biotype == "protein_coding")
PC_exon_pos <- PC_exon_pos[,c("chromosome_name", "exon_chrom_start", "exon_chrom_end")]
PC_exon_pos <- PC_exon_pos[complete.cases(PC_exon_pos),]
colnames(PC_exon_pos) <- c("chromosome", "start", "end")
# add splice AD sites to exon coordinates 
PC_exon_pos$start <- PC_exon_pos$start - 2
PC_exon_pos$end <- PC_exon_pos$end + 2

# subset NC exon POS
NC_exon_pos <- subset(exon_pos, exon_pos$gene_biotype != "protein_coding")
NC_exon_pos <- NC_exon_pos[,c("chromosome_name", "exon_chrom_start", "exon_chrom_end")]
NC_exon_pos <- NC_exon_pos[complete.cases(NC_exon_pos),]
colnames(NC_exon_pos) <- c("chromosome", "start", "end")
# add splice AD sites to exon coordinates 
NC_exon_pos$start <- NC_exon_pos$start - 2
NC_exon_pos$end <- NC_exon_pos$end + 2

# subsset UTR POS
UTR5 <- UTR_pos[,c("chromosome_name", "5_utr_start", "5_utr_end")]
UTR5 <- UTR5[complete.cases(UTR5),]
colnames(UTR5) <- c("chromosome", "start", "end")
UTR3 <- UTR_pos[,c("chromosome_name", "3_utr_start", "3_utr_end")]
UTR3 <- UTR3[complete.cases(UTR3),]
colnames(UTR3) <- c("chromosome", "start", "end")
UTR_pos <- rbind(UTR5, UTR3)

### subset intron POS
intron_pos <- ddply(exon_pos, "ensembl_transcript_id", get_intron_POS)
intron_pos$length <- intron_pos$end - intron_pos$start
long_intron <- subset(intron_pos, intron_pos$length >= 1000000)
intron_pos <- subset(intron_pos, intron_pos$length < 1000000)
intron_pos <- intron_pos[,c("chromosome", "start", "end")]
intron_pos <- intron_pos[complete.cases(intron_pos),]

# export data
fwrite(PC_exon_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/PC_exon_POS_short.csv")
fwrite(NC_exon_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/NC_exon_POS_short.csv")
fwrite(intron_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/intron_POS_short.csv")
fwrite(UTR_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/UTR_POS_short.csv")
fwrite(promoter_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/promoter_POS_short.csv")
fwrite(promoter_flank_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/promoter_flanking_POS_short.csv")
fwrite(enhancer_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/enhancer_POS_short.csv")
fwrite(TF_binding_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/TF_binding_POS_short.csv")
fwrite(open_chromatin_pos, "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/open_chromatin_POS_short.csv")


####
