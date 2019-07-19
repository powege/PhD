rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

### SET ARGS

# chromo <- args[1]
# ref.file <- args[2]
# annotated.file <- args[3]
# ranked.file <- args[4]
# out.file <- args[5]
chromo <- 19
ref.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Mouse_REF_sm_Ensembl_GRCm38_v94_chr19.txt"
annotated.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation.csv"
ranked.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
out.file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_annotation_CpG_proportion_mouse_chr19.csv"


### IMPORT

# import REF 
ref <- fread(ref.file)
# import annotation
ann <- fread(annotated.file)
ranked.ann<- fread(ranked.file)


### FORMAT

# identify REF start and end POS 
start.POS <- which(ref$REF != "N")[1]
end.POS <- which(ref$REF != "N")[length(which(ref$REF != "N"))]
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
# remove start < start.POS | end > end.POS
ann <- ann[which(ann$start >= start.POS & ann$end <= end.POS)]
# ann <- subset(ann, ann$start >= start.POS & ann$end <= end.POS)

### Get CpG proportion for each sequence

## OLD and SLOW 
# ann$CpG_proportion <- CpG_prop_vectorise(from = ann$start, to = ann$end)
# ann$CpG_proportion <- NA
# for(i in 1:nrow(ann)){
#   string <- paste(
#     ref$REF [ c( which(ref$POS == ann$start[i]) : which(ref$POS == ann$end[i]) ) ], 
#     collapse = "")
#   CpG_count <- str_count(string, "CG")
#   total.dinuc <- nchar(string)/2
#   ann$CpG_proportion[i] <- CpG_count/total.dinuc
#   print(i)
# }

# collapse vector to string
ref_string <- paste0(ref$REF, collapse="")

# get sequence for each annotation
ann$sequence <- sapply(1:nrow(ann), function(i) {
    substr(ref_string, ann$start[i], ann$end[i])
  })
# get CpG proportion for each annotation
ann$CpG_proportion <- sapply(1:nrow(ann), function(i) {
  str_count(ann$sequence[i], "CG")/(nchar(ann$sequence[i])/2)
})
# get annotation length
ann$length <- nchar(ann$sequence)
# ann$length <- (ann$end - ann$start) +1
# get proportion of "N"
ann$N_proportion <- str_count(ann$sequence, "N")


### EXPORT
out <- ann[,c("category", "chromosome", "CpG_proportion", "N_proportion", "length")]
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
