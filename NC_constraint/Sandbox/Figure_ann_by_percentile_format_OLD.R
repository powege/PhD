rm(list=ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} 

ind <- as.integer(args[1])
annotation.file <- args[2]
ranked.file <- args[3]
out.file <- args[4]
species <- args[5]
# ind <- 1
# annotation.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_annotation.csv"
# unannotated.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRCh38_GENCODE_RegBuild_unannotated.csv"
# out.file <- "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/Annotation_territory_by_percentile_650_50.csv"
# species <- "human"

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


# IMPORT
ann <- fread(annotation.file)
int_un <- fread(ranked.file)

# subset Intron and Unannotated from ranked.file
int_un <- subset(int_un, int_un$category == "Intron" | int_un$category == "Unannotated")

# remove Intron from annotation.file
ann <- subset(ann, ann$category != "Intron")

# rbind annotated and unannotated
dt <- rbind(ann, int_un)
rm(ann, int_un)

### FORMAT

# list all annotations
all_annotations <- c(
  "Exon - CDS",
  "Exon - UTR", 
  "Exon - non-coding", 
  "Promoter", 
  "Enhancer",
  "Open chromatin",
  "TF binding",
  "Promoter flanking", 
  "Intron",
  "Unannotated"
)
annotation <- all_annotations[ind]

# subset annotation
dt_sub <- subset(dt, dt$category == annotation)

# get vector of all annotation POS_ID 
out <- list()
if (species == "human") { chr <- c(1:22) }
if (species == "mouse") { chr <- c(1:19) }
for(i in chr){
  sub <- subset(dt_sub, dt_sub$chromosome == chr[i])
  all.POS <- unlist(seq2(from = sub$start, to = sub$end))
  all.CHR <- rep(chr[i], length(all.POS))
  out[[i]] <- paste0(all.CHR, "_", all.POS)
  print(i)
}
out <- unique(unlist(out))


# get proportion of annotation POS in each percentile
prop <- rep(NA, 100)
percentile <- c(1:100)
for (i in percentile){
  p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/", species, "_POS_ID_percentile_", i, "_650_50.csv"), header = F)
  prop[i] <- sum(p$V1 %chin% out)/length(p$V1)
  print(i)
}

df_out <- data.frame(Percentile = percentile, 
                 Fraction = prop,
                 Annotation = rep(annotation, 100))

### EXPORT
fwrite(df_out, out.file, append = T, col.names = F)


# ### FOR INTRONS
# 
# if (annotation == "Intron"){
#   
#   dt_sub <- subset(dt, dt$category == annotation)
#   
#   if (species == "human") { 
#     chr <- c(1:22) 
#   for(i in chr){
#   sub <- subset(dt_sub, dt_sub$chromosome == i)
#   all.POS <- unlist(seq2(from = sub$start, to = sub$end))
#   all.CHR <- rep(i, length(all.POS))
#   out <- paste0(all.CHR, "_", all.POS)
#   out <- as.data.table(unique(out)) # ensure unique
#   fwrite(
#     out, 
#     "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/human_intron_POS_ID_long.csv",
#     append = T,
#     col.names = F
#   )
#   print(i)
#   }
#   rm(out)
#   }
#   
#   if (species == "mouse") { 
#     chr <- c(1:19) 
#     for(i in chr){
#       sub <- subset(dt_sub, dt_sub$chromosome == i)
#       all.POS <- unlist(seq2(from = sub$start, to = sub$end))
#       all.CHR <- rep(i, length(all.POS))
#       out <- paste0(all.CHR, "_", all.POS)
#       out <- as.data.table(unique(out)) # ensure unique
#       fwrite(
#         out, 
#         "/well/lindgren/George/Data/NC_constraint/Ensembl_annotation_POS/mouse_intron_POS_ID_long.csv",
#         append = T,
#         col.names = F
#       )
#       print(i)
#     }
#     rm(out)
#   }
#   
# }

#####

### STACK OVERFLOW

# # how to improve performance?
# set.seed(17)
# vec1 <- paste0(sample(1:10, 10000000, replace = T), "_", sample(1:1000000000, 10000000))
# vec2 <- paste0(sample(1:10, 1000000, replace = T), "_", sample(1:1000000000, 1000000))
# 
# system.time({ prop1 <- table(vec2 %in% vec1)[[2]]/length(vec2) })
# system.time({ prop2 <- length(intersect(vec2, vec1)) / length(vec2) })
# system.time({ prop3 <- sum(vec2 %chin% vec1)/length(vec2) })

# I have two very large vectors that I need to concatenate with a delimiter to form unique IDs. For example: 
#   
#   set.seed(1)
# 
# vec1 <- sample(1:10, 10000000, replace = T)
# vec2 <- sample(1:1000000000, 10000000))
# 
# I am currently using paste0():
#   
#   system.time({    
#     
#     uniq_id <- paste0(vec1, "_", vec2)
#     
#   })
# 
# However, due to the size of vec1 and vec2 this is quite slow. Is there an alternate method with greater performance? 
