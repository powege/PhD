### SCRIPT that formats GENCODE and Ensembl Regulatory Build annotation
### Outputs cols: category; chromosome; start; end

rm(list=ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are three argument: if not, return an error
if (length(args)!=3) {
  stop("Exactly three arguments must be supplied", call.=FALSE)
} 

# set args variables
reg_build.gff <- args[1] # .gff file
gencode.gtf <- args[2] # .gtf file
out_file <- args[3] # .csv file
# reg_build.gff <- "~/Dropbox/PhD/Data/Ensembl/Annotation/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff"
# gencode.gtf <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Homo_sapiens.GRCh38.94.gtf"
# out_file <- "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/Human_GENCODE_RegBuild_annotation.csv"

### FUNCTIONS
# sub <- subset(exon_all, exon_all$transcript_id == "ENST00000450305")
### FUNCTION that infers intron POS from exon POS
get_intron_POS <- function(sub){
  if (nrow(sub) == 1){
    out <- data.frame(chromosome = sub$chromosome,
                      start = NA,
                      end = NA)
  }
  if (nrow(sub) > 1){
    sub <- sub[order(sub$exon_number),]
    start <- sub$end[1:(nrow(sub))-1] + 1
    end <- sub$start[2:(nrow(sub))] - 1
    chrom <- rep(sub$chromosome[1], length(start))
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


### IMPORT annotation raw data
reg_build <-  fread(reg_build.gff)
gene <- fread(gencode.gtf)
# motif <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/dummy/Homo_sapiens.GRCh38.motif_features.gff.dummy")


### FORMAT gene annotation
colnames(gene) <- c("chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
# pull out transcript_biotype and transcript_id from attribute
# gene$gene_biotype <- str_match(string=gene$attribute, pattern="gene_biotype\\s+\"([^\"]+)\"")[,2]
gene$transcript_biotype <- str_match(string=gene$attribute, pattern="transcript_biotype\\s+\"([^\"]+)\"")[,2]
gene$transcript_id <- str_match(string=gene$attribute, pattern="transcript_id\\s+\"([^\"]+)\"")[,2]
gene$exon_number <- str_match(string=gene$attribute, pattern="exon_number\\s+\"([^\"]+)\"")[,2]
# subset CDS and UTR
cds <- subset(gene, gene$type == "CDS")
utr <- subset(gene, gene$type == "five_prime_utr" | gene$type == "three_prime_utr")
# identify transcript_id in CDS and UTR 
PCtranscript <- unique(c(cds$transcript_id, utr$transcript_id))
# subset all exons
exon_all <- subset(gene, gene$type == "exon")
# subset all exon with non-PC transcript
exon_nonPC <- exon_all[which(!exon_all$transcript_id %in% PCtranscript),]
# identify introns from exons
intron_all <- ddply(exon_all, "transcript_id", get_intron_POS)


### FORMAT Regulatory Build
colnames(reg_build) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
# split attribute column
reg_build <-
  reg_build %>% separate(attribute, c("ID", "bound_end","bound_start", "description", "feature_type"), sep = ";")
# remove tags 
reg_build$ID <- gsub("ID=", "", reg_build$ID)
reg_build$bound_end <- gsub("bound_end=", "", reg_build$bound_end)
reg_build$bound_start <- gsub("bound_start=", "", reg_build$bound_start)
reg_build$description <- gsub("description=", "", reg_build$description)
reg_build$feature_type <- gsub("feature_type=", "", reg_build$feature_type)


### FORMAT annotated categories

# subset promotor POS
promoter_pos <- subset(reg_build, reg_build$feature_type == "Promoter")
promoter_pos$category <- "Promoter"
promoter_pos <- promoter_pos[,c("category","chromosome", "start", "end")]
promoter_pos <- promoter_pos[complete.cases(promoter_pos),]
# subset promotor flanking POS
promoter_flank_pos <- subset(reg_build, reg_build$feature_type == "Promoter Flanking Region")
promoter_flank_pos$category <- "Promoter flanking"
promoter_flank_pos <- promoter_flank_pos[,c("category", "chromosome", "start", "end")]
promoter_flank_pos <- promoter_flank_pos[complete.cases(promoter_flank_pos),]
# subset enhancer POS
enhancer_pos <- subset(reg_build, reg_build$feature_type == "Enhancer")
enhancer_pos$category <- "Enhancer"
enhancer_pos <- enhancer_pos[,c("category", "chromosome", "start", "end")]
enhancer_pos <- enhancer_pos[complete.cases(enhancer_pos),]
# subset TF binding site POS
TF_binding_pos <- subset(reg_build, reg_build$feature_type == "TF binding site" | reg_build$feature_type == "CTCF Binding Site")
TF_binding_pos$category <- "TF binding"
TF_binding_pos <- TF_binding_pos[,c("category", "chromosome", "start", "end")]
TF_binding_pos <- TF_binding_pos[complete.cases(TF_binding_pos),]
# subset open chromatin POS
open_chromatin_pos <- subset(reg_build, reg_build$feature_type == "Open chromatin")
open_chromatin_pos$category <- "Open chromatin"
open_chromatin_pos <- open_chromatin_pos[,c("category", "chromosome", "start", "end")]
open_chromatin_pos <- open_chromatin_pos[complete.cases(open_chromatin_pos),]
# subset exon CDS POS
cds$category <- "Exon - CDS"
exon_CDS_pos <- cds[,c("category", "chromosome", "start", "end")]
exon_CDS_pos <- exon_CDS_pos[complete.cases(exon_CDS_pos),]
# subset exon UTR POS
utr$category <- "Exon - UTR"
exon_UTR_pos <- utr[,c("category", "chromosome", "start", "end")]
exon_UTR_pos <- exon_UTR_pos[complete.cases(exon_UTR_pos),]
# subset exon non-coding POS
exon_nonPC$category <- "Exon - non-coding"
exon_NC_pos <- exon_nonPC[,c("category", "chromosome", "start", "end")]
exon_NC_pos <- exon_NC_pos[complete.cases(exon_NC_pos),]
# subset intron POS
intron_all$category <- "Intron"
intron_pos <- intron_all[,c("category", "chromosome", "start", "end")]
intron_pos <- intron_pos[complete.cases(intron_pos),]
# rbind
annotation <- rbind(exon_CDS_pos,
                    exon_UTR_pos,
                    exon_NC_pos,
                    promoter_pos,
                    promoter_flank_pos,
                    enhancer_pos,
                    TF_binding_pos,
                    open_chromatin_pos,
                    intron_pos)
# remove duplicates
annotation <- unique(annotation)
# annotation <- annotation[!duplicated(annotation),]

### EXPORT 
fwrite(annotation, out_file)


### RESOURCES 
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md



