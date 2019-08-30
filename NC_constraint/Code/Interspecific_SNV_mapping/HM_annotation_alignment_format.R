rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


#######################################################################
### Q1: What fraction of the genome aligns between human and mouse? ###
#######################################################################

# Edit so as not to include centromeres and telomeres!!!!

## X% of human genome aligns to mouse
# import all human seqences that align to mouse
h_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_short.txt")
# set human chr lengths (GRCh38) from 
h_chr_len <- c(248956422,	242193529,198295559,190214555,181538259,170805979,159345973,145138636,
               138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
               83257441,80373285,58617616,64444167,46709983,50818468)
# ,156040895)

# get alignment fraction for each autosome
h_align_frac_chr <- rep(NA, 22)
h_chr <- c(1:22)
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

## import annotations
# human
h_ann_r <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation_ranked.csv")
h_ann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv")
# mouse
m_ann_r <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation_ranked.csv")
m_ann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv")

## import alignment
h_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/H_HtoM_alignment_short.txt")
m_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/M_HtoM_alignment_short.txt")


# combine unannotated from ranked
h_ann <- subset(h_ann, h_ann$category != "Intron")
h_ann_r <- subset(h_ann_r, h_ann_r$category == "Intron" | h_ann_r$category == "Unannotated")
h_ann <- rbind(h_ann, h_ann_r)
m_ann <- subset(m_ann, m_ann$category != "Intron")
m_ann_r <- subset(m_ann_r, m_ann_r$category == "Intron" | m_ann_r$category == "Unannotated")
m_ann <- rbind(m_ann, m_ann_r)

# ensure start <= end 
tmp1 <- subset(h_ann, h_ann$end < h_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(h_ann, h_ann$end >= h_ann$start)
h_ann <- rbind(tmp1, tmp2)
tmp1 <- subset(m_ann, m_ann$end < m_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(m_ann, m_ann$end >= m_ann$start)
m_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)

h_dt_out <- data.table()
h_chr <- c(1:22)
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





h_dt_out <- data.table()
h_chr <- c(1:22)
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



# categorise all mouse to human alignment by mouse annotation

# M_chr; M_start; M_end; H_chr; H_start; H_end; M_anotation; H_annotation

# H_chr; H_start; H_end; H_annotation

# for each mouse annotation
  # identify all positions that align to human
    # H_chr; H_start; H_end; H_annotation


# Human

# annotation; chr; start; end; 
# total_align; Exon_CDS_align; Exon_UTR_align;
# Exon_noncoding_align; Promoter_align; Promoter_flanking_align;
# Enhancer_align; Open_chromatin_align; TF_binding_align;
# Intron_align; Unannotated_align



