rm(list = ls())
graphics.off()

library(data.table)

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT
h_align <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_short.csv")
h_ann <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv")


# set human chr lengths (GRCh38) from 
h_chr_len <- c(248956422,	242193529,198295559,190214555,181538259,170805979,159345973,145138636,
               138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
               83257441,80373285,58617616,64444167,46709983,50818468,156040895)
h_chr <- c(1:22, "X")
m_chr <- c(1:19, "X")

### FORMAT

colnames(h_ann) <- c("category", "chromosome", "start", "end")
colnames(h_align) <- c("H_chr", "H_start", "H_end", "H_ann", "M_chr", "M_start", "M_end", "M_ann")
h_align <- subset(h_align, h_align$M_chr %in% m_chr)
h_align <- h_align[,c("H_chr", "H_start", "H_end")]

# ensure start <= end 
tmp1 <- subset(h_ann, h_ann$end < h_ann$start)
colnames(tmp1) <- c("category", "chromosome", "end", "start")
tmp2 <- subset(h_ann, h_ann$end >= h_ann$start)
h_ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)
tmp1 <- subset(h_align, h_align$H_end < h_align$H_start)
colnames(tmp1) <- c("H_chr",  "H_end", "H_start") 
tmp2 <- subset(h_align, h_align$H_end >= h_align$H_start)
h_align <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)


### What fraction of each genomic annotation aligns between human and mouse? 
# annotations can overlap (ie one base can have multipe annotations, including of the same type)

h_dt_out <- data.table()
for(i in 1:length(h_chr)){
  df1_sub <- subset(h_ann, h_ann$chromosome == h_chr[i])
  ann_vec <- df1_sub$category
  df1_sub <- df1_sub[,c('start', 'end')]
  df2_sub <- subset(h_align, h_align$H_chr == h_chr[i])
  df2_int <- unique(unlist(seq2(from = df2_sub$H_start, to = df2_sub$H_end)))
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

h_dt_out$length <- abs(h_dt_out$end - h_dt_out$start) + 1
h_dt_out$align_frac <- h_dt_out$MATCH/h_dt_out$length

### EXPORT
fwrite(h_dt_out, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_annotation_alignment_by_sequence.csv")
# h_dt_out <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation.csv")


#######################################################################
### Q1: What fraction of the genome aligns between human and mouse? ###
#######################################################################

### ACCOUNT FOR OVERLAP!
### ACCOUNT FOR UNALIGNED MASKED REGIONS! 

# # get alignment fraction for each autosome
# h_align_frac_chr <- rep(NA, 22)
# h_chr <- c(1:22)
# for (i in 1:length(h_chr)){
#   tmp <- subset(h_align, h_align$H_chr == h_chr[i])
#   h_align_frac_chr[i] <- sum((tmp$H_end + 1) - tmp$H_start) / h_chr_len[i]
#   print(i)
# }
# # get total autosome alignment
# h_align_frac <- sum(h_align_frac_chr * h_chr_len[1:22]) / sum(h_chr_len[1:22])


