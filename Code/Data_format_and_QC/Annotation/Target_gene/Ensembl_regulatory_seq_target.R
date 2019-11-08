### SCRIPT that returns closest genes to a given sequence
# format regulatory input
# format TSS domains
# all canonical TSS with overlap
# max one canonical TSS per promoter prioritised by proximal overlap

rm(list=ls())
graphics.off()

library(data.table)
library( matrixStats )
library(plyr)

### FUNCTIONS

dominion <- function(tss, distal_extension){
  
# subset by strand 
strand1 <- subset(tss, tss$strand == 1)
strand2 <- subset(tss, tss$strand == -1)

# get start and end pos 
strand1$proximal_start <- strand1$transcription_start_site - 5000
strand1$proximal_end <- strand1$transcription_start_site + 1000
strand1$distal1_start <- strand1$transcription_start_site - (distal_extension + 5000)
strand1$distal1_end <- strand1$transcription_start_site - 5001
strand1$distal2_start <- strand1$transcription_start_site + 1001
strand1$distal2_end <- strand1$transcription_start_site + (distal_extension + 1000)

strand2$proximal_start <- strand2$transcription_start_site + 5000
strand2$proximal_end <- strand2$transcription_start_site - 1000
strand2$distal1_start <- strand2$transcription_start_site + (distal_extension + 5000)
strand2$distal1_end <- strand2$transcription_start_site + 5001
strand2$distal2_start <- strand2$transcription_start_site - 1001
strand2$distal2_end <- strand2$transcription_start_site - (distal_extension + 1000)

# rbind
tss <- rbind(strand1, strand2)
rm(strand1, strand2)

# melt 
id_vars <- c("chromosome","external_gene_name","ensembl_gene_id","ensembl_transcript_id","transcription_start_site",
             "transcript_biotype","strand")
dt_start <- tss[,c(id_vars,"proximal_start", "distal1_start","distal2_start"), with = F]
dt_start <- melt(dt_start, id.vars = id_vars)
colnames(dt_start) <- c(id_vars, "domain", "domain_start")
dt_start$domain <- as.character(dt_start$domain)
dt_start$domain[dt_start$domain == "proximal_start"] <- "proximal"
dt_start$domain[dt_start$domain == "distal1_start"] <- "distal1"
dt_start$domain[dt_start$domain == "distal2_start"] <- "distal2"

dt_end <- tss[,c(id_vars,"proximal_end", "distal1_end","distal2_end"), with = F]
dt_end <- melt(dt_end, id.vars = id_vars)
colnames(dt_end) <- c(id_vars, "domain", "domain_end")
dt_end$domain <- as.character(dt_end$domain)
dt_end$domain[dt_end$domain == "proximal_end"] <- "proximal"
dt_end$domain[dt_end$domain == "distal1_end"] <- "distal1"
dt_end$domain[dt_end$domain == "distal2_end"] <- "distal2"

tss <- dt_start[dt_end, on = c(id_vars, "domain")]
tss <- unique(tss)
rm(dt_start, dt_end)

return(tss)
}

all_canonical <- function(reg, tss, chr, can_ID){
  
  # subset chromosome
  reg_sub <- subset(reg, reg$chromosome == chr)
  tss_sub <- subset(tss, tss$chromosome == chr)
  
  # merge by overlap
  setkey(reg_sub, chromosome, category_start, category_end)
  setkey(tss_sub, chromosome, domain_start, domain_end)
  ans <- unique( foverlaps( tss_sub, reg_sub ))
  
  # calculate total overlap
  ans[, overlap_start := rowMaxs( as.matrix(.SD), na.rm = TRUE ), .SDcols = c("category_start", "domain_start")]
  ans[, overlap_end   := rowMins( as.matrix(.SD), na.rm = TRUE ), .SDcols = c("category_end", "domain_end")]
  ans[, overlap_size  := overlap_end - overlap_start + 1 ]
  
  # complete cases
  ans <- ans[complete.cases(ans),]
  
  # calculate min distance to TSS from cat start and end
  ans$TSS_dist <- apply( data.table( start_dist = abs( (ans$transcription_start_site + 1) - ans$category_start ),
                                     end_dist = abs( (ans$transcription_start_site + 1) - ans$category_end) ), 
                         1,
                         min
  )
  
  # overlap fractions
  ans$category_length <- (ans$category_end + 1) - ans$category_start
  ans$domain_length <- (ans$domain_end + 1) - ans$domain_start
  ans$f_category_overlap <- ans$overlap_size/ans$category_length
  ans$f_domain_overlap <- ans$overlap_size/ans$domain_length
  
  # subset variables
  ans <- ans[,c("chromosome", "category_ID", "category", "category_start", "category_end", "external_gene_name",
                "ensembl_gene_id", "ensembl_transcript_id", "transcription_start_site",
                "domain", "TSS_dist", "f_category_overlap", "f_domain_overlap")]
  ans <- ans[order(category, category_start, category_end),] # order 
  ans <- ans[ensembl_transcript_id %in% can_ID] # subset canonical
  
  return(ans)
}

# FUNCTION that returns row with shortest TSS_dist
# sub <- out_p[category_ID == "D1"]
closest_one <- function(sub){ sub[which(sub$TSS_dist == min(sub$TSS_dist)),] }

# FUNCTION that prioritises proximal overlap
prioritise_proximal <- function(sub){ 
  if ("proximal" %in% sub$domain) {
    out <- sub[which(sub$domain == "proximal"),]
  } else ( out <- sub )
  return(out)
}

# FUNCTION that prioritises one promoter per transcript by TSS distance
# sub <- out_p[out_p$ensembl_transcript_id == "ENST00000247977",]
one_canonical <- function(sub){
  if ("proximal" %in% sub$domain) {
    out <- sub[which(sub$domain == "proximal" & sub$TSS_dist == min(sub$TSS_dist)),]
  } else ( out <- sub[which(sub$TSS_dist == min(sub$TSS_dist)),] )
}

# FUNCTION that runs script
alakazam <- function(tss_file, reg_file, can_file, n_chr){
  
  ### IMPORT
  
  # TSS
  tss <- fread(tss_file)
  # annotation sequences
  reg <- fread(reg_file)
  # canonical transcripts
  can <- fread(can_file)
  
  
  ### SET VARS 
  
  long_order <- c("Exon - CDS",
                  "Exon - UTR",
                  "Exon - other",
                  "Promoter",
                  "Enhancer - proximal",
                  "Enhancer - distal",
                  "CTCF binding",
                  "Miscellaneous",
                  "Intron - proximal",
                  "Intron - distal",
                  "Unannotated")
  code_order <- c("A","B","C","D","E","F","G","H","I","J","K")
  can_ID <- unique(can$ensembl_transcript_id) # unique canonical IDs 
  rm(can)
  
  
  ### FORMAT
  
  ## format regulatory input

  # colnames
  colnames(reg) <- c("chromosome", "category_start", "category_end", "category", "category_strand")
  # convert miscellaneous
  reg$category[reg$category == "TF binding" | reg$category == "Open chromatin"] <- "Miscellaneous"
  # subset regulatory elements
  reg <- reg[category %in% c("Promoter", "Enhancer - distal", "Enhancer - proximal", "CTCF binding")]
  # ensure start >= end
  reg[ category_start > category_end, `:=`( category_start = category_end, category_end = category_start)]
  # add annotation ID
  reg$category_code <- mapvalues(reg$category, from=long_order, # number code annotations (plyr)
                                 to=code_order)
  tmp_list <- list()
  for (cat in 1:length(unique(reg$category))){
    tmp_sub <- subset(reg, reg$category == unique(reg$category)[cat])
    tmp_sub$category_ID <- paste0(tmp_sub$category_code, 1:length(tmp_sub$category))
    tmp_list[[cat]] <- tmp_sub
  }
  reg <- do.call("rbind", tmp_list)
  reg$category_code <- NULL
  
  ## format TSS domains
  
  # colnames
  colnames(tss) <- c("chromosome",
                     "external_gene_name",
                     "ensembl_gene_id",
                     "ensembl_transcript_id",
                     "transcription_start_site",
                     "transcript_biotype",
                     "strand")
  # proximal and distal domains
  tss <- dominion(tss, distal_extension = 300000)
  # ensure start >= end
  tss[ domain_start > domain_end, `:=`( domain_start = domain_end, domain_end = domain_start)]
  # overlap for all canonical TSS
  chr_list <- list()
  for(chr in 1:n_chr){
    chr_list[[chr]] <- all_canonical(reg = reg, tss = tss, chr = chr, can_ID = can_ID)
    print(chr)
  }
  output_all <- do.call("rbind", chr_list)
  # one canonical transcript per promoter
  out_p <- output_all[category == "Promoter"]
  out_p <- ddply(out_p, "category_ID", closest_one)
  out_p <- ddply(out_p, "category_ID", prioritise_proximal)
  out_p <- ddply(out_p, "ensembl_transcript_id", one_canonical)
  output <- rbind(out_p, output_all[category %in% c("Enhancer - distal", "Enhancer - proximal", "CTCF binding")])
  
  return(output)
}

### SCRIPT

# run for human
h_out <- alakazam(tss_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_TSS.csv",
         reg_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv",
         can_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_human_canPC_pos.csv",
         n_chr = 22)
fwrite(h_out, "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_target_canTSS.csv")

# run for mouse
m_out <- alakazam(tss_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_TSS.csv",
         reg_file = "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation.csv",
         can_file = "~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_mouse_canPC_pos.csv",
         n_chr = 19)
fwrite(m_out, "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRC38_GENCODE_RegBuild_annotation_target_canTSS.csv")






#####

# ### STACK OVERFLOW
# 
# ### QUESTION
# dt1 <- data.table(ID1 = c("A", "B", "C", "D", "E"),
#                   start1 = c(100, 1, 210, 300, 400),
#                   end1 = c(200, 90, 240, 380, 500))
# dt2 <- data.table(ID2 = c("a1", "a2", "a3", "a4", "a5", "a6"),
#                   start2 = c(10, 150, 300, 310, 350, 400),
#                   end2 = c(50, 100, 250, 280, 390, 450))
# output <- data.table(ID1 = c("A", "B", "D", "D", "D", "E"),
#                      start1 = c(100, 1, 300, 300, 300, 400),
#                      end1 = c(200, 90, 380, 380, 380, 500),
#                      ID2 = c("a2", "a1", "a3", "a4", "a5", "a6"),
#                      start2 = c(150, 10, 300, 310, 350, 400),
#                      end2 = c(100, 50, 250, 280, 390, 450))
# 
# ID1_list <- list() # set output lists 
# ID2_list <- list()
# for (i in 1:nrow(dt1)){
#   vec1 <- seq(from = dt1$start1[i], to = dt1$end1[i])
#   ID1_vec <- rep(dt1$ID1, each = nrow(dt2)) # set output vectors
#   ID2_vec <- rep(NA, nrow(dt2))
#   for (j in 1:nrow(dt2)){
#     vec2 <- seq(from = dt2$start[j], to = dt2$end[j])
#     if (length(intersect(vec2, vec1)) > 0){
#       ID2_vec[j] <- dt2$ID2[j]
#     }
#   }
#   ID1_list[[i]] <- ID1_vec
#   ID2_list[[i]] <- ID2_vec
# }
# output2 <- data.table(ID1 = unlist(ID1_list),
#                       ID2 = unlist(ID2_list))
# output2 <- output2[complete.cases(output2),]
# output2 <- merge(dt1, unique(output2))
# output2 <- merge(output2, dt2, by = "ID2")
# 
# ### ANSWER
# library(data.table)
# 
# # tmp1.1 <- subset(dt1, dt1$end1 >= dt1$start1)
# # tmp1.2 <- subset(dt1, dt1$end1 < dt1$start1)
# # colnames(tmp1.2) <- c("ID1", "end1", "start1")
# # dt1 <- rbind(tmp1.1, tmp1.2)
# # rm(tmp1.1, tmp1.2)
# # tmp2.1 <- subset(dt2, dt2$end2 >= dt2$start2)
# # tmp2.2 <- subset(dt2, dt2$end2 < dt2$start2)
# # colnames(tmp2.2) <- c("ID2", "end2", "start2")
# # dt2 <- rbind(tmp2.1, tmp2.2)
# # rm(tmp2.1, tmp2.2)
# 
# #in foverlaps(), start should always be before end..
# #so switch dt2's values where this is not the case
# dt2[ start2 > end2, `:=`( start2 = end2, end2 = start2)]
# 
# setkey(dt1, start1, end1)
# setkey(dt2, start2, end2)
# output3 <- foverlaps( dt2, dt1 )
# 
# for (i in 1:nrow(output3)){
#   tmp_vec1 <- seq(from = output3$start1[i], to = output3$end1[i]) 
#   tmp_vec2 <- seq(from = output3$start2[i], to = output3$end2[i]) 
#   output3$overlap[i] <- length(intersect(tmp_vec1, tmp_vec2))
# }
# 
# ans <- foverlaps( dt2, dt1 )
# library( matrixStats )
# ans[, overlap_start := rowMaxs( as.matrix(.SD), na.rm = TRUE ), .SDcols = c("start1", "start2")]
# ans[, overlap_end   := rowMins( as.matrix(.SD), na.rm = TRUE ), .SDcols = c("end1", "end2")]
# ans[, overlap_size  := overlap_end - overlap_start + 1 ]
# 
# #####
# 
# df1 <- data.table(CAT = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#                   START = c(1, 11, 21, 1, 21, 41, 1, 11, 21),
#                   END = c(10, 20, 30, 20, 40, 60, 10, 20, 30)
# )
# df2 <- data.table(CAT = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#                   START = c(1, 11, 21, 31, 41, 51, 1, 11, 21),
#                   END = c(5, 17, 23, 38, 48, 54, 9, 17, 26)
# )
# foverlaps(df1[, rn := .I], setkey(df2, CAT, START, END))[
#   , ovl := (pmin(END, i.END) - pmax(START, i.START) + 1)][
#     , .(MATCH = sum(ovl)), by = .(rn)][
#       is.na(MATCH), MATCH := 0][]
