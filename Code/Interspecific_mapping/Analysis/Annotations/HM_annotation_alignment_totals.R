rm(list=ls())
graphics.off()

library(data.table)
library(plyr)

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there are three argument: if not, return an error
# if (length(args)<1) {
#   stop("More that one argument must be supplied", call.=FALSE)
# } 

### set args variables
# ann_align_prefix  <- args[1]
# H_ann_prefix <- args[2]
# out_file <- args[3]

ann_align_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_short.csv"
H_ann_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation.csv"
out_file <- "~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_annotation_alignment_totals.csv"

# ann_align_prefix  <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_short_Hchr"
# H_ann_prefix <- "/well/lindgren/George/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_chr"
# out_file <- "/well/lindgren/George/Data/Ensembl/Alignment/Formatted/HM_alignment_annotation_summary.csv"

### FUNCTIONS
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# set vars
h_chr <- c(1:22)
m_chr <- c(1:19, "X")
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
dt_out <- data.table()

### IMPORT
dt_Hann_all <- fread(H_ann_file)
dt_align_all <- fread(ann_align_file)

### FORMAT

colnames(dt_Hann_all) <- c("category", "chromosome", "start", "end")
dt_Hann_all$category[dt_Hann_all$category == "TF binding" | dt_Hann_all$category == "Open chromatin"] <- "Miscellaneous"
dt_Hann_all$category <- mapvalues(dt_Hann_all$category, from=long_order, # number code annotations (plyr)
                              to=code_order)
# dt_Hann_all$len <- abs(dt_Hann_all$end - dt_Hann_all$start) + 1 # length

colnames(dt_align_all) <- c("H_chromosome", "H_start", "H_end", "H_category", "M_chromosome", "M_start", "M_end", "M_category")
dt_align_all <- subset(dt_align_all, dt_align_all$M_chromosome %in% m_chr) # subset Mouse chr
# dt_align_all$len <- abs(dt_align_all$H_start - dt_align_all$H_end) + 1 # sequence length

### FOR LOOP
for (chr in 1:length(h_chr)){

# subset chr
dt_Hann_sub <- subset(dt_Hann_all, dt_Hann_all$chromosome == h_chr[chr])
dt_align_sub <- subset(dt_align_all, dt_align_all$H_chromosome == h_chr[chr])

# set output objects
H_chromosome <- rep(h_chr[chr], length(code_order)*length(code_order))
H_annotation <- list()
M_annotation <- list()
H_ann_total <- list()
H_ann_align_total <- list()
HM_ann_align_total <- list()
tmp_vec <- rep(NA, length(code_order))

for (i in 1:length(code_order)){
  
  # subset annotation
  sub <- subset(dt_Hann_sub, dt_Hann_sub$category == code_order[i])
  # total annotation bp (no overlap)
  H_ann_total[[i]] <- rep(length(unique(as.vector(unlist(seq2(from = sub$start, to = sub$end))))), length(code_order))
  
  # subset annotation
  sub1 <- subset(dt_align_sub, dt_align_sub$H_category %like% code_order[i])
  # total annotation bp that align (no overlap)
  H_ann_align_total[[i]] <- rep(length(unique(as.vector(unlist(seq2(from = sub1$H_start, to = sub1$H_end))))), length(code_order))
  
  for (j in 1:length(code_order)){
    
    # subset annotation
    sub2 <- subset(sub1, sub1$M_category %like% code_order[j])
    # total annotation bp that align to the same annotation in mouse (no overlap)
    tmp_vec[j] <- length(unique(as.vector(unlist(seq2(from = sub2$H_start, to = sub2$H_end)))))
    
    print(c(i,j))
  }
  
  HM_ann_align_total[[i]] <- tmp_vec
  H_annotation[[i]] <- rep(code_order[i], length(code_order))
  M_annotation[[i]] <- code_order
}  
  

dt_tmp <- data.table(H_chromosome = H_chromosome,
                     H_annotation = unlist(H_annotation),
                     H_ann_total = unlist(H_ann_total),
                     H_ann_align_total = unlist(H_ann_align_total),
                     M_annotation = unlist(M_annotation),
                     HM_ann_align_total = unlist(HM_ann_align_total)
                     )

dt_tmp$H_annotation <- mapvalues(dt_tmp$H_annotation, from=code_order, to=long_order)
dt_tmp$M_annotation <- mapvalues(dt_tmp$M_annotation, from=code_order, to=long_order)

dt_out <-  rbind(dt_out, dt_tmp)

print(h_chr[chr])
}

dt_out$ind <- paste0(dt_out$H_annotation, dt_out$M_annotation)

# sub <- subset(dt_out, dt_out$ind == "PromoterExon - CDS")
format_fun <- function(sub){
  data.table(H_annotation = sub$H_annotation[1],
             H_ann_total = sum(sub$H_ann_total),
             H_ann_align_total = sum(sub$H_ann_align_total),
             M_annotation = sub$M_annotation[1],
             HM_ann_align_total = sum(sub$HM_ann_align_total)
  )
}

dt_plot <- ddply(dt_out, "ind", format_fun)
dt_plot$ind <- NULL
dt_plot$H_ann_align_frac <- dt_plot$H_ann_align_total/dt_plot$H_ann_total
dt_plot$HM_ann_align_frac <- dt_plot$HM_ann_align_total/dt_plot$H_ann_align_total 

### EXPORT
fwrite(dt_plot, out_file)


#####
