rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


#######################################################################
### Q1: What fraction of the genome aligns between human and mouse? ###
#######################################################################

## X% of human genome aligns to mouse
# import all human seqences that align to mouse
h_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_short.txt")
# set human chr lengths (GRCh38) from 
h_chr_len <- c(248956422,	242193529,198295559,190214555,181538259,170805979,159345973,145138636,
               138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
               83257441,80373285,58617616,64444167,46709983,50818468,156040895)
# get alignment fraction for each autosome
h_align_frac_chr <- rep(NA, 23)
h_chr <- c(1:22, "X")
for (i in 1:length(h_chr)){
  tmp <- subset(h_align, h_align$V1 == h_chr[i])
  h_align_frac_chr[i] <- sum((tmp$V3 + 1) - tmp$V2) / h_chr_len[i] 
}
# get total autosome alignment
h_align_frac <- sum(h_align_frac_chr * h_chr_len) / sum(h_chr_len)


################################################################################
### Q2: What fraction of each genomic annotation aligns between human and mouse? ###
################################################################################
  
### For each human annotation, what percentage aligns to mouse?

# import annotations
h_ann_r <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv")
h_ann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv")
# import alignment
h_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_short.txt")

# combine unannotated from ranked
h_ann <- subset(h_ann, h_ann$category != "Intron")
h_ann_r <- subset(h_ann_r, h_ann_r$category == "Intron" | h_ann_r$category == "Unannotated")
h_ann <- rbind(h_ann, h_ann_r)

# ensure start <= end 
tmp1 <- subset(h_ann, h_ann$end < h_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(h_ann, h_ann$end >= h_ann$start)
h_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)

h_dt_out <- data.table()
h_chr <- c(1:22, "X")
for(i in 1:length(h_chr)){
  df1_sub <- subset(h_ann, h_ann$chromosome == h_chr[i])
  ann_vec <- df1_sub$category
  df1_sub <- df1_sub[,c('start', 'end')]
  df2_sub <- subset(h_align, h_align$V1 == h_chr[i])
  df2_int <- unique(unlist(seq2(from = df2_sub$V2, to = df2_sub$V3)))
  df1_sub[, ind := .I] # add uniqe index to data.table
  dt2_sub <- as.data.table(df2_int, key = 'vec1') # convert to data.table
  colnames(dt2_sub) <- 'vec1'
  dt2_sub[, vec2 := df2_int] # dublicate column
  setkey(df1_sub) # sets keys // order data by all columns
  # Fast overlap join:
  ans1 = foverlaps(dt2_sub, df1_sub, by.x = c('vec1', 'vec2'), by.y = c('start', 'end'),
                   type = "within", nomatch = 0L)
  
  counts <- ans1[, .N, keyby = ind] # count by ind
  # merge to inital data
  df1_sub[, MATCH := counts[df1_sub, on = .(ind), x.N]]
  
  setorder(df1_sub, ind) # reorder by ind to get inital order
  df1_sub[, ind := NULL] # deletes ind colum
  df1_sub[is.na(MATCH), MATCH := 0L] # NAs is 0 count
  df1_sub$chromosome <- h_chr[i]
  df1_sub$category <- ann_vec
  h_dt_out <- rbind(h_dt_out, df1_sub)
  
  print(h_chr[i])
}
h_dt_out$length <- (h_dt_out$end + 1) - h_dt_out$start
h_dt_out$align_frac <- h_dt_out$MATCH/h_dt_out$length
fwrite(h_dt_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation.csv")

# totals
cat_vec <- unique(h_dt_out$category)
len_vec <- rep(NA, length(cat_vec)) 
match_vec <- rep(NA, length(cat_vec)) 
align_frac_vec <- rep(NA, length(cat_vec)) 
for (i in 1:length(cat_vec)){
  tmp_sub <- subset(h_dt_out, h_dt_out$category == cat_vec[i])
  len_vec[i] <- sum(as.numeric(tmp_sub$length))
  match_vec[i] <- sum(as.numeric(tmp_sub$MATCH))
  align_frac_vec[i] <- match_vec[i]/len_vec[i]
}
h_dt_total_out <- data.table(annotation = cat_vec,
                           ann_bp_total = len_vec,
                           ann_bp_aligned = match_vec,
                           ann_aligned_frac = align_frac_vec)
fwrite(h_dt_total_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation_totals.csv")

### What is the annotation composition of the alignment? 

###################################################################################
### Q3: What fraction of human SNVs align to mouse (all gnomAD, ClinVar, GWAS)? ###
###################################################################################

### ClinVar
CV <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_formatted.vcf")
CV <- subset(CV, CV$CAT != "")
CV <- CV[,c("CHR", "POS", "REF", "CAT")]
CV <- CV[!duplicated(CV),] # removes SNVs with multiple ALT
  
CV_align <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_human_to_mouse_alignment.txt")
CV_align <- CV[CV_align, on = c("CHR", "POS", "REF")]
CV_align <- CV_align[,c("CHR", "POS", "REF", "CAT")]
CV_align <- CV_align[!duplicated(CV_align),]

CV <- as.data.frame(table(CV$CAT))
colnames(CV) <- c("annotation", "n_SNV")
CV_align <- as.data.frame(table(CV_align$CAT))
colnames(CV_align) <- c("annotation", "n_SNV_align")

h_dt_total_out <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation_totals.csv")
CV_out <- merge(CV, CV_align)
CV_out$n_SNV_align_frac <- CV_out$n_SNV_align/CV_out$n_SNV
CV_out <- merge(CV_out, h_dt_total_out[,c("annotation", "ann_aligned_frac")])

### OUTPUT
fwrite(CV_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_syntenic.csv")


### GWAS
GWAS <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_formatted.txt")
GWAS <- subset(GWAS, GWAS$CAT != "")
GWAS <- GWAS[,c("CHR", "POS", "CAT")]
GWAS <- GWAS[!duplicated(GWAS),] # removes SNVs with multiple ALT

GWAS_align <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_human_to_mouse_alignment.txt")
GWAS_align <- GWAS[GWAS_align, on = c("CHR", "POS")]
GWAS_align <- GWAS_align[,c("CHR", "POS", "REF", "CAT")]
GWAS_align <- GWAS_align[!duplicated(GWAS_align),]

GWAS <- as.data.frame(table(GWAS$CAT))
colnames(GWAS) <- c("annotation", "n_SNV")
GWAS_align <- as.data.frame(table(GWAS_align$CAT))
colnames(GWAS_align) <- c("annotation", "n_SNV_align")

h_dt_total_out <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation_totals.csv")
GWAS_out <- merge(GWAS, GWAS_align)
GWAS_out$n_SNV_align_frac <- GWAS_out$n_SNV_align/GWAS_out$n_SNV
GWAS_out <- merge(GWAS_out, h_dt_total_out[,c("annotation", "ann_aligned_frac")])

### OUTPUT
fwrite(GWAS_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_syntenic.csv")



#########

GWAS_out$CAT <- "GWAS"
CV_out$CAT <- "ClinVar"
x <- rbind(CV_out, GWAS_out)


# ## X% of mouse genome aligns to human
# # import all mouse seqences that align to human
# m_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/M_HtoM_alignment_short.txt")
# # set human chr lengths (GRCh38) from 
# m_chr_len <- c(195471971, 182113224, 160039680, 156508116,	151834684,	
#                149736546, 145441459,	129401213, 124595110,	130694993,	
#                122082543, 120129022,	120421639,	124902244,	104043685,	
#                98207768,	94987271,	90702639,	61431566, 171031299)
# # get alignment fraction for each autosome
# m_align_frac_chr <- rep(NA, 20)
# m_chr <- c(1:19, "X")
# for (i in 1:length(m_chr)){
#   tmp <- subset(m_align, m_align$V1 == m_chr[i])
#   m_align_frac_chr[i] <- sum((tmp$V3 + 1) - tmp$V2) / m_chr_len[i] 
# }
# # get total autosome alignment
# m_align_frac <- sum(m_align_frac_chr * m_chr_len) / sum(m_chr_len)
# 
# align_frac_chr <- data.table(species = c(rep("human", 24), rep("mouse", 21)),
#                              chr = c(1:22, "X", "total", 1:19, "X", "total"),
#                              length = c(h_chr_len, sum(h_chr_len), m_chr_len, sum(m_chr_len)),
#                              fraction = c(h_align_frac_chr, h_align_frac, m_align_frac_chr, m_align_frac))
# 
# fwrite(align_frac_chr, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_chr.csv")

# ## Mouse
# 
# m_dt_out <- data.table()
# m_chr <- c(1:19, "X")
# for(i in 1:length(m_chr)){
#   df1_sub <- subset(m_ann, m_ann$chromosome == m_chr[i])
#   ann_vec <- df1_sub$category
#   df1_sub <- df1_sub[,c('start', 'end')]
#   df2_sub <- subset(m_align, m_align$V1 == m_chr[i])
#   df2_int <- unique(unlist(seq2(from = df2_sub$V2, to = df2_sub$V3)))
#   df1_sub[, ind := .I] # add uniqe index to data.table
#   dt2_sub <- as.data.table(df2_int, key = 'vec1') # convert to data.table
#   colnames(dt2_sub) <- 'vec1'
#   dt2_sub[, vec2 := df2_int] # dublicate column
#   setkey(df1_sub) # sets keys // order data by all columns
#   # Fast overlap join:
#   ans1 = foverlaps(dt2_sub, df1_sub, by.x = c('vec1', 'vec2'), by.y = c('start', 'end'),
#                    type = "within", nomatch = 0L)
#   
#   counts <- ans1[, .N, keyby = ind] # count by ind
#   # merge to inital data
#   df1_sub[, MATCH := counts[df1_sub, on = .(ind), x.N]]
#   
#   setorder(df1_sub, ind) # reorder by ind to get inital order
#   df1_sub[, ind := NULL] # deletes ind colum
#   df1_sub[is.na(MATCH), MATCH := 0L] # NAs is 0 count
#   df1_sub$chromosome <- m_chr[i]
#   df1_sub$category <- ann_vec
#   m_dt_out <- rbind(m_dt_out, df1_sub)
#   
#   print(m_chr[i])
# }
# m_dt_out$length <- (m_dt_out$end + 1) - m_dt_out$start
# m_dt_out$align_frac <- m_dt_out$MATCH/m_dt_out$length
# 
# # totals
# cat_vec <- unique(m_dt_out$category)
# len_vec <- rep(NA, length(cat_vec)) 
# match_vec <- rep(NA, length(cat_vec)) 
# align_frac_vec <- rep(NA, length(cat_vec)) 
# for (i in 1:length(cat_vec)){
#   tmp_sub <- subset(m_dt_out, m_dt_out$category == cat_vec[i])
#   len_vec[i] <- sum(tmp_sub$length)
#   match_vec[i] <- sum(tmp_sub$MATCH)
#   align_frac_vec[i] <- match_vec[i]/len_vec[i]
# }
# m_dt_total_out <- data.table(annotation = cat_vec,
#                              ann_bp_total = len_vec,
#                              ann_bp_aligned = match_vec,
#                              ann_aligned_frac = align_frac_vec)
# fwrite(m_dt_total_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_mouse_annotation.csv")



### PLOT 

GWAS_plot <- GWAS_out
GWAS_plot$n_SNV_align_exp <- GWAS_plot$n_SNV * GWAS_plot$ann_aligned_frac
GWAS_plot$n_SNV_frac <- 1

# GWAS_plot <- GWAS_plot[,c("annotation", "n_SNV", "n_SNV_align", "n_SNV_align_exp")]
# colnames(GWAS_plot) <- c("annotation", "Human", "Mouse", "Mouse expected")
# GWAS_plot <- melt(GWAS_plot, by = "annotation")

GWAS_plot <- GWAS_plot[,c("annotation", "n_SNV_frac", "n_SNV_align_frac", "ann_aligned_frac")]
colnames(GWAS_plot) <- c("annotation", "Human", "Mouse", "Syntenic SNVs")
GWAS_plot <- melt(GWAS_plot, by = "annotation")

p_gwas <- ggplot(data=GWAS_plot, aes(x=annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("GWAS") +
  xlab("Huaman annotation") +
  coord_flip()
p_gwas

CV_plot <- CV_out
CV_plot$n_SNV_align_exp <- CV_plot$n_SNV * CV_plot$ann_aligned_frac
CV_plot$n_SNV_frac <- 1

# CV_plot <- CV_plot[,c("annotation", "n_SNV", "n_SNV_align", "n_SNV_align_exp")]
# colnames(CV_plot) <- c("annotation", "Human", "Mouse", "Mouse expected")
# CV_plot <- melt(CV_plot, by = "annotation")

CV_plot <- CV_plot[,c("annotation", "n_SNV_frac", "n_SNV_align_frac", "ann_aligned_frac")]
colnames(CV_plot) <- c("annotation", "Human", "Mouse", "Syntenic SNVs")
CV_plot <- melt(CV_plot, by = "annotation")

p_cv <- ggplot(data=CV_plot, aes(x=annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("ClinVar") +
  xlab("Huaman annotation") +
  coord_flip()
p_cv

pout <- grid.arrange(p_cv, p_gwas, nrow = 1, widths = c(5, 5))


GWAS_plot <- GWAS_out
GWAS_plot$n_SNV_align_exp <- GWAS_plot$n_SNV * GWAS_plot$ann_aligned_frac
GWAS_plot$n_SNV_frac <- 1
CV_plot <- CV_out
CV_plot$n_SNV_align_exp <- CV_plot$n_SNV * CV_plot$ann_aligned_frac
CV_plot$n_SNV_frac <- 1

# ##################
# ### STACK OVERFLOW
# ##################
# 
# # ann
# df1 <- data.table(CAT = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#                   START = c(1, 11, 21, 1, 21, 41, 1, 11, 21),
#                   END = c(10, 20, 30, 20, 40, 60, 10, 20, 30)
# )
# # align
# df2 <- data.table(CAT = c(rep("A", 3), rep("B", 3), rep("C", 3)),
#                   START = c(1, 11, 21, 31, 41, 51, 1, 11, 21),
#                   END = c(5, 17, 23, 38, 48, 54, 9, 17, 26)
# )
# 
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
# df1$MATCH <- NA
# for (i in 1:nrow(df1)){
#   df2_sub <- subset(df2, df2$CAT == df1$CAT[i])
#   df2_int <- unlist(seq2(from = df2_sub$START, to = df2_sub$END))
#   df1_int <- seq(df1$START[i], df1$END[i])
#   df1$MATCH[i] <- length(na.omit(match(df1_int, df2_int)))
#   print(i)
# }
# 
# dt_out <- data.table()
# for(i in 1:length(unique(df1$CAT))){
#   df1_sub <- subset(df1, df1$CAT == unique(df1$CAT)[i])
#   df1_sub <- df1_sub[,c('START', 'END')]
#   df2_sub <- subset(df2, df1$CAT == unique(df1$CAT)[i])
#   df2_int <- unlist(seq2(from = df2_sub$START, to = df2_sub$END))
#   df1_sub[, ind := .I] # add uniqe index to data.table
#   dt2_sub <- as.data.table(df2_int, key = 'vec1') # convert to data.table
#   colnames(dt2_sub) <- 'vec1'
#   dt2_sub[, vec2 := df2_int] # dublicate column
#   setkey(df1_sub) # sets keys // order data by all columns
#   # Fast overlap join:
#   ans1 = foverlaps(dt2_sub, df1_sub, by.x = c('vec1', 'vec2'), by.y = c('START', 'END'),
#                    type = "within", nomatch = 0L)
#   
#   counts <- ans1[, .N, keyby = ind] # count by ind
#   # merge to inital data
#   df1_sub[, MATCH := counts[df1_sub, on = .(ind), x.N]]
#   
#   setorder(df1_sub, ind) # reorder by ind to get inital order
#   df1_sub[, ind := NULL] # deletes ind colum
#   df1_sub[is.na(MATCH), MATCH := 0L] # NAs is 0 count
#   df1_sub$CAT <- unique(df1$CAT)[i]
#   dt_out <- rbind(dt_out, df1_sub)
#   
#   print(unique(df1$CAT)[i])
# }

# colnames(m_align) <- c("chromosome", "start", "end")
# m_ann$chromosome <- as.character(m_ann$chromosome)
# foverlaps(m_ann[, rn := .I], setkey(m_align, chromosome, start, end))[
#   , ovl := (pmin(end, i.end) - pmax(start, i.start) + 1)][
#     , .(MATCH = sum(ovl)), by = .(rn)][
#       is.na(MATCH), MATCH := 0][]

