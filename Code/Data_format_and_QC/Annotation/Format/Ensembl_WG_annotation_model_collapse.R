### SCRIPT that formats Ensembl gene and multicell Regulatory Build annotation
### merges seqences with overlapping annotations of same category
### inferes Intron and Unannotated regions by chromosome
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

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
reg_build.gff <- args[1] # .gff file
gencode.gtf <- args[2] # .gtf file
out_file <- args[3] # .csv file
chr <- args[4]
species <- args[5]
# reg_build.gff <- "~/Dropbox/PhD/Data/Ensembl/Annotation/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff"
# gencode.gtf <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Homo_sapiens.GRCh38.94.gtf"
# out_file <- "~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/Human_GENCODE_RegBuild_annotation.csv"
# chr <- as.integer(17)
# species <- "homo_sapien"


### FUNCTIONS

### FUNCTION that collapses overlapping sequences
# sub <- subset(promoter_pos, promoter_pos$chromosome == "3")
collapse_overlap <- function(sub){
  
  # ensure start <= end
  tmp <- subset(sub, sub$start > sub$end)
  colnames(tmp) <- c("category", "chromosome", "end", "start")
  dt <- rbind(subset(sub, sub$start <= sub$end), tmp)

  dt2 <- rbind(dt,dt)
  dt2[1:(.N/2), ID := 1]
  dt2[(.N/2 +1):.N, ID := 2]
  setkey(dt2, ID, start, end)
  
  squished <- dt2[,.(START_DT = start, 
                     END_DT = end, 
                     indx = c(0, cumsum(as.numeric(lead(start)) > cummax(as.numeric(end)))[-.N])),
                  keyby=ID
                  ][,.(start=min(START_DT), 
                       end = max(END_DT)),
                    by=c("ID","indx")
                    ]
  squished <- squished[, c("start", "end")]
  squished <- unique(squished)
  squished$category <- sub$category[1]
  squished$chromosome <- sub$chromosome[1]
  
  return(squished[,c("category", "chromosome", "start", "end")])
}


# sub <- subset(exon_all, exon_all$transcript_id == "ENST00000623180")
### FUNCTION that infers intron - distil POS from exon POS
get_intron_proximal_POS <- function(sub){
  if (nrow(sub) == 1){
    out <- data.frame(chromosome = sub$chromosome,
                      start = NA,
                      end = NA)
  }
  if (nrow(sub) > 1){
    sub <- sub[order(sub$end),] # order by exon end POS!!!
    start <- c((sub$end[1:(nrow(sub))-1] + 1), (sub$start[2:(nrow(sub))] - 10))
    end <- c((sub$end[1:(nrow(sub))-1] + 10), (sub$start[2:(nrow(sub))] - 1))
    chrom <- rep(sub$chromosome[1], length(start))
    out <- data.frame(chromosome = chrom,
                      start = start,
                      end = end)
  }
  return(out)
}

# sub <- subset(exon_all, exon_all$transcript_id == "ENST00000450305")
### FUNCTION that infers intron POS from exon POS
get_intron_POS <- function(sub){
  if (nrow(sub) == 1){
    out <- data.frame(chromosome = sub$chromosome,
                      start = NA,
                      end = NA)
  }
  if (nrow(sub) > 1){
    # sub <- sub[order(sub$exon_number),] # order by exon number
    sub <- sub[order(sub$end),] # order by exon end POS!!!
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

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# GRCh38 chromosome lengths
H_chr <- c(1:22)
H_chr_length <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
                  145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                  101991189,  90338345,  83257441,  80373285,  58617616,  64444167,  46709983,
                  50818468)
M_chr <- c(1:19)
M_chr_length <- c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 
                  129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 
                  104043685, 98207768, 94987271, 90702639, 61431566)
if (species == "homo_sapien"){ 
  chr_length <- H_chr_length
  chromosomes <- H_chr }
if (species == "mus_musculus"){ 
  chr_length <- M_chr_length 
  chromosomes <- M_chr }

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
promoter_flank_pos$category <- "Enhancer - proximal"
promoter_flank_pos <- promoter_flank_pos[,c("category", "chromosome", "start", "end")]
promoter_flank_pos <- promoter_flank_pos[complete.cases(promoter_flank_pos),]
# subset enhancer POS
enhancer_pos <- subset(reg_build, reg_build$feature_type == "Enhancer")
enhancer_pos$category <- "Enhancer - distal"
enhancer_pos <- enhancer_pos[,c("category", "chromosome", "start", "end")]
enhancer_pos <- enhancer_pos[complete.cases(enhancer_pos),]
# subset CTCF binding site POS
CTCF_binding_pos <- subset(reg_build, reg_build$feature_type == "CTCF Binding Site")
CTCF_binding_pos$category <- "CTCF binding"
CTCF_binding_pos <- CTCF_binding_pos[,c("category", "chromosome", "start", "end")]
CTCF_binding_pos <- CTCF_binding_pos[complete.cases(CTCF_binding_pos),]
# subset TF binding site POS
TF_binding_pos <- subset(reg_build, reg_build$feature_type == "TF binding site")
TF_binding_pos$category <- "Miscellaneous"
TF_binding_pos <- TF_binding_pos[,c("category", "chromosome", "start", "end")]
TF_binding_pos <- TF_binding_pos[complete.cases(TF_binding_pos),]
# subset open chromatin POS
open_chromatin_pos <- subset(reg_build, reg_build$feature_type == "Open chromatin")
open_chromatin_pos$category <- "Miscellaneous"
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
exon_nonPC$category <- "Exon - other"
exon_NC_pos <- exon_nonPC[,c("category", "chromosome", "start", "end")]
exon_NC_pos <- exon_NC_pos[complete.cases(exon_NC_pos),]

# rbind
annotation <- rbind(exon_CDS_pos,
                    exon_UTR_pos,
                    exon_NC_pos,
                    promoter_pos,
                    promoter_flank_pos,
                    enhancer_pos,
                    CTCF_binding_pos,
                    TF_binding_pos,
                    open_chromatin_pos)
# remove duplicates
annotation <- unique(annotation)
annotation <- annotation[complete.cases(annotation),]
# annotation <- annotation[!duplicated(annotation),]


### FORMAT introns 

# subset all exons
exon_all <- subset(exon_all, exon_all$chromosome == chr)
# identify intron - proximal from exons
intron_prox_pos <- ddply(exon_all, "transcript_id", get_intron_proximal_POS)
intron_prox_pos <- intron_prox_pos[complete.cases(intron_prox_pos),]
intron_prox_pos$category <- "Intron - proximal"
intron_prox_pos$transcript_id <- NULL
# subset annotation by chr
annotation <- subset(annotation, annotation$chromosome == chr)
# rbind 
annotation <- rbind(annotation, intron_prox_pos)
# identify intron - distil from exons
intron_dis_all <- ddply(exon_all, "transcript_id", get_intron_POS)
intron_dis_all <- intron_dis_all[complete.cases(intron_dis_all),]
# exclude regions already annotated
tmp_ann <- annotation[,c("start", "end")]
colnames(tmp_ann) <- c("from", "to")
tmp_int <- intron_dis_all[,c("start", "end")]
colnames(tmp_int) <- c("from", "to")
# vector of intron POS not in other annotations
intron_pos <- setdiff(x = unlist(seq2(from = tmp_int$from, to = tmp_int$to)), 
                      y = unlist(seq2(from = tmp_ann$from, to = tmp_ann$to)))
# intron_pos <- setdiff(x = 1:chr_length[chr], y = unlist(pmap(tmp, seq)))
# format vector to dataframe
intron_pos <- sort(intron_pos, decreasing = F)
intron_pos <- t(sapply(split(intron_pos, findInterval(intron_pos, intron_pos[which(c(1, diff(intron_pos)) > 1)])), range))
intron_pos <- as.data.table(intron_pos)
colnames(intron_pos) <- c("start", "end")
intron_pos$chromosome <- chr
intron_pos$category <- "Intron - distal"
intron_pos <- unique(intron_pos)
# rbind 
annotation <- rbind(annotation, intron_pos)


### FORMAT unannotated

tmp_ann <- annotation[,c("start", "end")]
colnames(tmp_ann) <- c("from", "to")
# vector of intron POS not in other annotations
unan_pos <- setdiff(x = 1:chr_length[as.integer(chr)], 
                      y = unlist(seq2(from = tmp_ann$from, to = tmp_ann$to)))
# unan_pos <- setdiff(x = 1:chr_length[chr], y = unlist(pmap(tmp, seq)))
# format vector to dataframe
unan_pos <- sort(unan_pos, decreasing = F)
unan_pos <- t(sapply(split(unan_pos, findInterval(unan_pos, unan_pos[which(c(1, diff(unan_pos)) > 1)])), range))
unan_pos <- as.data.table(unan_pos)
colnames(unan_pos) <- c("start", "end")
unan_pos$chromosome <- chr
unan_pos$category <- "Unannotated"
unan_pos <- unique(unan_pos)
# rbind 
annotation <- rbind(annotation, unan_pos)


### Colapse overlapping sequences for autosomes
tmp_list <- list()
for (i in 1:length(unique(annotation$category))){
    tmp_list[[i]] <- collapse_overlap(subset(annotation, annotation$category == unique(annotation$category)[i]))
}
annotation <- do.call("rbind", tmp_list)


### EXPORT 
fwrite(annotation, out_file, col.names = F)


### RESOURCES 
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


#####

### STACK OVERFLOW

# dt <- data.table(start = c(1, 11, 20, 44, 3, 46, 52),
#                  end = c(8, 15, 13, 50, 9, 46, 52))
# 
# out <- data.table(start = c(1, 11, 44),
#                   end = c(9, 20, 50))
# 
# tmp1 <- subset(dt, dt$start > dt$end)
# colnames(tmp1) <- c("end", "start")
# dt <- rbind(subset(dt, dt$start <= dt$end),
#             tmp1)
# dt <- dt[order(start),]
# 
# dt2 <- rbind(dt,dt)
# dt2[1:(.N/2), ID := 1]
# dt2[(.N/2 +1):.N, ID := 2]
# setkey(dt2, ID, start, end)
# 
# squished <- dt2[,.(START_DT = start, 
#                   END_DT = end, 
#                   indx = c(0, cumsum(as.numeric(lead(start)) > cummax(as.numeric(end)))[-.N])),
#                keyby=ID
#                ][,.(start=min(START_DT), 
#                     end = max(END_DT)),
#                  by=c("ID","indx")
#                  ]
# squished <- squished[, c("start", "end")]
# squished <- unique(squished)
# 
# start_col = c("2018-01-01","2018-03-01","2018-03-10","2018-03-20","2018-04-10","2018-05-01","2018-05-05","2018-05-10","2018-07-07")
# end_col = c("2018-01-21","2018-03-21","2018-03-31","2018-04-09","2018-04-30","2018-05-21","2018-05-26","2018-05-30","2018-07-14")
# 
# d <- data.table(start_col = as.Date(start_col), end_col = as.Date(end_col))
# d2<- rbind(d,d)
# d2[1:(.N/2), ID := 1]
# d2[(.N/2 +1):.N, ID := 2]
# d2[17,end_col := as.Date('2018-12-01')]
# 
# setkey(d2, ID, start_col, end_col)
# 
# squished <- d2[,.(START_DT = start_col, 
#                   END_DT = end_col, 
#                   indx = c(0, cumsum(as.numeric(lead(start_col)) > cummax(as.numeric(end_col)))[-.N])),
#                keyby=ID
#                ][,.(start=min(START_DT), 
#                     end = max(END_DT)),
#                  by=c("ID","indx")
#                  ]



