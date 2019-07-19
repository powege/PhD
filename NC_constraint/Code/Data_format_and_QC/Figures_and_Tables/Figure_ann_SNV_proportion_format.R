rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

### SET ARGS

chromo <- args[1]
vcf.file <- args[2]
annotated.file <- args[3]
ranked.file <- args[4]
out.file <- args[5]
# chromo <- 19
# vcf.file <- "~/Dropbox/PhD/Data/MGP/vcf_QCed_VEP/MGP_v5_snps_QCed_hom_chr19.vcf"
# annotated.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv"
# ranked.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
# out.file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_annotation_SNV_proportion_mouse_chr19.csv"


### IMPORT

# import SNV file
vcf <- fread(vcf.file)
# import annotation
ann <- fread(annotated.file)
ranked.ann<- fread(ranked.file)


### FORMAT

# subset ann != "Intron"
ann <- subset(ann, ann$category != "Intron")
# subset ranked Intron and Unannotated
ranked.ann <- subset(ranked.ann, ranked.ann$category == "Intron" | ranked.ann$category == "Unannotated")
# rbind files
ann <- rbind(ann, ranked.ann)
rm(ranked.ann)
# subset chromo
ann <- subset(ann, ann$chromosome == chromo)
# ensure start <= end
tmp1 <- subset(ann, ann$start <= ann$end)
tmp2 <- subset(ann, ann$start > ann$end)
colnames(tmp2) <- c("category", "chromosome", "end", "start")
ann <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)

### Get fraction of SNVs for each annotation

# COUNT <- rep(NA, nrow(ann)) 
# for (i in 1:nrow(ann)){
#   tmp_vec <- seq(from = ann$start[i], to = ann$end[i])
#   COUNT[i] <- table(tmp_vec %in% vcf$V2)[2]
#   print(i)
# }
# ann$n_SNV <- COUNT

ann[, ind := .I] # add uniqe index to data.table
tmp_ann <- ann[,c("start", "end")]
tmp_ann[, ind := .I] # add uniqe index to data.table
tmp_dt <- as.data.table(vcf$V2) # convert to data.table
tmp_dt[, V2 := V1] # dublicate column
setkey(tmp_ann) # sets keys // order data by all columns
# Fast overlap join:
ans1 = foverlaps(tmp_dt, tmp_ann, by.x = c('V1', 'V2'), by.y = c('start', 'end'),
                 type = "within", nomatch = 0L)

counts <- ans1[, .N, keyby = ind] # count by ind
# merge to inital data
ann[, n_SNV := counts[ann, on = .(ind), x.N]]

setorder(ann, ind) # reorder by ind to get inital order
ann[, ind := NULL] # deletes ind colum
ann[is.na(n_SNV), n_SNV := 0L] # NAs is 0 count

# get length of each annotation in kb
ann$length_kb <- ((ann$end - ann$start) + 1 )/1000

# get adjusted n_SNV
ann$n_SNV_kb <- ann$n_SNV / ann$length_kb

### EXPORT
out <- ann[,c("category", "chromosome", "n_SNV", "length_kb", "n_SNV_kb")]
fwrite(out, out.file, append = T)


#####

# ### FUNCTIONS
# 
# # FUNCTION that calculates the proportion of CpG dinucleotides in sequence given start and end coordinates
# CpG_prop <- function(from, to){
#   string <- paste(
#     ref$REF [ c( which(ref$POS == from) : which(ref$POS == to) ) ], 
#     collapse = "")
#   CpG_count <- str_count(string, "CG")
#   total.dinuc <- nchar(string)/2
#   CpG_prop <- CpG_count/total.dinuc
#   return(CpG_prop)
# }
# 
# # Function that vectorises CpG_prop
# CpG_prop_vectorise <- Vectorize(CpG_prop, vectorize.args = c("from", "to"))
# 
# 
# # ### SPEED TEST
# 
# dummy <- ann[38001:38100,]
# 
# # vectorise
# system.time({
#   test1 <- CpG_prop_vectorise(from = dummy$start, to = dummy$end)
#   }) # 50.838
# 
# # for loop
# system.time({
# test2 <- rep(NA, nrow(dummy))
# for(i in 1:nrow(dummy)){
#   string <- paste(
#     ref$REF [ c( which(ref$POS == dummy$start[i]) : which(ref$POS == dummy$end[i]) ) ],
#     collapse = "")
#   CpG_count <- str_count(string, "CG")
#   total.dinuc <- nchar(string)/2
#   test2[i] <- CpG_count/total.dinuc
#   print(i)
# }
# }) # 50.585

# ### STACK OVERFLOW
# 
# set.seed(1)
# df1 <- data.frame(POS = 1:10000000,
#                   REF = sample(c("A", "T", "G", "C"), 10000000, replace = T))
# 
# df2 <- data.frame(start = sample(1:5000000, 10, replace = T),
#                      end = sample(5000001:10000000, 10, replace = T))
# 
# ### For loop
# 
# system.time( {
# for(i in 1:nrow(df2)){
#   string <- paste(
#     df1$REF [ c( which(df1$POS == df2$start[i]) : which(df1$POS == df2$end[i]) ) ],
#     collapse = "")
#     CpG_count <- str_count(string, "CG")
#     total.dinuc <- nchar(string)/2
#     df2$CpG_proportion[i] <- CpG_count/total.dinuc
#   print(i)
# }
# })
# 
# # FUNCTION that calculates the proportion of CpG dinucleotides in sequence given start and end coordinates
# mongoose <- function(from, to){
#   string <- paste(
#     df1$REF [ c( which(df1$POS == from) : which(df1$POS == to) ) ],
#     collapse = "")
#   return(string)
# }
# 
# # Function that vectorises CpG_prop
# mongoose_vec <- Vectorize(mongoose, vectorize.args = c("from", "to"))
# 
# system.time({
#   sequences <- mongoose_vec(from = df2$start, to = df2$end)
#   })
# 
# 
# #
# ref2 <- paste0(df1$REF, collapse="")
# system.time({
#   sequences2 <- sapply(1:nrow(df2), function(i) {
#     substr(ref2, df2$start[i], df2$end[i])
#   })
# })
# CpG2 <- sapply(1:nrow(df2), function(i) {
#   str_count(sequences2[i], "CG")/(nchar(sequences2[i])/2)
# })

### STACK OVERFLOW

# set.seed(1)

# df1 <- data.table(
#   START = seq(1, 10000000, 10),
#   END = seq(10, 10000000, 10)
#   )
#                                           
# vec1 <- sample(1:100000, 10000)

# COUNT <- rep(NA, nrow(df1)) 
# for (i in 1:nrow(df1)){
#   vec2 <- seq(from = df1$START[i], to = df1$END[i])
#   COUNT[i] <- table(vec2 %in% vec1)[2]
#   print(i)
# }
# df1$COUNT <- COUNT

# test1 <- df1[, count := sum(between(vec1, START, END)), by = seq_len(nrow(df1))]

# df1[, ind := .I] # add uniqe index to data.table
# dt2 <- as.data.table(vec1, key = 'vec1') # convert to data.table
# dt2[, vec2 := vec1] # dublicate column
# setkey(df1) # sets keys // order data by all columns
# # Fast overlap join:
# ans1 = foverlaps(dt2, df1, by.x = c('vec1', 'vec2'), by.y = c('START', 'END'),
#                  type = "within", nomatch = 0L)
# 
# counts <- ans1[, .N, keyby = ind] # count by ind
# # merge to inital data
# df1[, COUNT := counts[df1, on = .(ind), x.N]]
# 
# setorder(df1, ind) # reorder by ind to get inital order
# df1[, ind := NULL] # deletes ind colum
# df1[is.na(COUNT), COUNT := 0L] # NAs is 0 count


