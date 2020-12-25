rm(list=ls())
graphics.off()

library(data.table)

# INPUT:
# human annotation bed 
# mouse annotation bed
# human mouse alignment bed
# chromosome

# OUTPUT:
# bp counts for alignment and conservation

### SET VARS 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments
if (length(args)==0) {
  stop("Arguments must be supplied", call.=FALSE)
} 

# set args variables
h.ann.file <- args[1]
m.ann.file <- args[2]
align.file <- args[3]
chr <- as.character(args[4])
functions.file <- args[5]
out.all.file <- args[6]
out.specific.file <- args[7]

# h.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/h.sap_GRC38_v101_whole_genome_features_multicell.csv"
# m.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/m.mus_GRC38_v101_whole_genome_features_multicell.csv"
# align.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/hsap_grch38_v_mmus_grcm38_v101_alignment_chr15.bed"
# chr <- as.character(15)
# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Interspecific_mapping/v1/FUNCTIONS.R"
# out.all.file <- ""
# out.specific.file <- ""

### FUNCTIONS
source(functions.file)

### IMPORT 
human_ann <- fread(h.ann.file)
mouse_ann <- fread(m.ann.file)
align <- fread(align.file)

### FORMAT

# subset by human chromosome
human_ann <- human_ann[chromosome == chr]
align <- align[chromosome_human == chr]

# QC
human_ann <- human_ann[percentage_Nmask < 1]
mouse_ann <- mouse_ann[percentage_Nmask < 1]

### RUN FOR ALL ANNNOTATION

# for each human annotation in chromosome:
ann <- unique(human_ann$category)
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  h_ann_colapse <- collapse.overlap(human_ann[category == ann[i]][,c("chromosome", "start", "end")])
  n_total[i] <- sum(abs( (h_ann_colapse$end + 1) - h_ann_colapse$start))
  
  # orthologous sequences for human annotation to get total alignment
  ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                               "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = human_ann[category == ann[i]][,c("chromosome", "start", "end")])
  n_aligned[i] <- sum(abs( (ann_align$end_A + 1) - ann_align$start_A))

  # total of mouse annotation in orthologous sequences for conservation
  m_ann_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = mouse_ann[category == ann[i]][,c("chromosome", "start", "end")])
  n_conserved[i] <- sum(abs ((m_ann_conserved$end + 1) - m_ann_conserved$start))
  
  print(ann[i])
}

out_dt_all <- data.table( chromosome = chr,
                          annotation = ann,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved)
out_dt_all$p_aligned <- out_dt_all$n_aligned / out_dt_all$n_total
out_dt_all$p_conserved <- out_dt_all$n_conserved / out_dt_all$n_total
# rm(ann_align, h_ann_colapse, m_ann_conserved, ann, i, n_aligned, n_conserved, n_total)

### RUN FOR "Functional - all" "Non-CDS" and "Total"

human_ann_2_list <- list(functional = human_ann[category %in% c("Exon - CDS",
                                                           "Exon - UTR",
                                                           "Exon - other",
                                                           "Promoter",           
                                                           "Enhancer - proximal",
                                                           "Enhancer - distal",
                                                           "CTCF binding",
                                                           "Miscellaneous",      
                                                           "TAD boundry",
                                                           "Intron - proximal")],
                    non_cds =  human_ann[!category %in% c("Exon - CDS", "Intron - distal", "Unannotated")],
                    total = human_ann)

mouse_ann_2_list <- list(functional = mouse_ann[category %in% c("Exon - CDS",
                                                           "Exon - UTR",
                                                           "Exon - other",
                                                           "Promoter",           
                                                           "Enhancer - proximal",
                                                           "Enhancer - distal",
                                                           "CTCF binding",
                                                           "Miscellaneous",      
                                                           "TAD boundry",
                                                           "Intron - proximal")],
                    non_cds =  mouse_ann[!category %in% c("Exon - CDS", "Intron - distal", "Unannotated")],
                    total = mouse_ann)

ann <- c("Functional - all", "Functional - nonCDS", "Total")
n_total <- rep(NA, length(ann))
n_aligned <- rep(NA, length(ann))
n_conserved <- rep(NA, length(ann))
for(i in 1:length(ann)){
  
  human_ann_2 <- human_ann_2_list[[i]]
  mouse_ann_2 <- mouse_ann_2_list[[i]]
  
  h_ann_colapse <- collapse.overlap(human_ann_2[,c("chromosome", "start", "end")])
  n_total[i] <- sum(abs( (h_ann_colapse$end + 1) - h_ann_colapse$start))
  
  # orthologous sequences for human annotation to get total alignment
  ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                    "chromosome_mouse", "start_mouse", "end_mouse")],
                               bed = h_ann_colapse)
  n_aligned[i] <- sum(abs((ann_align$end_A + 1) - ann_align$start_A))
  sum(abs( (ann_align$end_B + 1) - ann_align$start_B))
  
  # total of mouse annotation in orthologous sequences for conservation
  m_ann_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                   bed2 = mouse_ann_2[,c("chromosome", "start", "end")])
  n_conserved[i] <- sum(abs( (m_ann_conserved$end + 1) - m_ann_conserved$start))
  
  print(ann[i])
}

out_dt_2 <- data.table( chromosome = chr,
                          annotation = ann,
                          n_total = n_total,
                          n_aligned = n_aligned,
                          n_conserved = n_conserved)
out_dt_2$p_aligned <- out_dt_2$n_aligned / out_dt_2$n_total
out_dt_2$p_conserved <- out_dt_2$n_conserved / out_dt_2$n_total

### RUN FOR TISSUE SPECIFIC ANNNOTATION

human_enhancer <- human_ann[category == "Enhancer - distal" | category == "Enhancer - proximal"]
mouse_enhancer <- mouse_ann[category == "Enhancer - distal" | category == "Enhancer - proximal"]
human_enhancer_specific <- list(heart = human_enhancer[activity_heart == "ACTIVE" | activity_heart == "POISED"],
                           kidney = human_enhancer[activity_kidney == "ACTIVE" | activity_kidney == "POISED"],
                           spleen = human_enhancer[activity_spleen == "ACTIVE" | activity_spleen == "POISED"])
mouse_enhancer_specific <- list(heart = mouse_enhancer[activity_heart == "ACTIVE" | activity_heart == "POISED"],
                           kidney = mouse_enhancer[activity_kidney == "ACTIVE" | activity_kidney == "POISED"],
                           spleen = mouse_enhancer[activity_spleen == "ACTIVE" | activity_spleen == "POISED"])

out_list <- list()
for (j in 1:length(human_enhancer_specific)){
  
  # for each human annotation in chromosome:
  human_enhancer_sub <- human_enhancer_specific[[j]]
  mouse_enhancer_sub <- mouse_enhancer_specific[[j]]
  
  ann <- unique(human_enhancer_sub$category)
  n_total <- rep(NA, length(ann))
  n_aligned <- rep(NA, length(ann))
  n_conserved <- rep(NA, length(ann))
  for(i in 1:length(ann)){
    
    # collapse annotation to get total bases
    h_enhancer_colapse <- collapse.overlap(human_enhancer_sub[category == ann[i]][,c("chromosome", "start", "end")])
    n_total[i] <- sum(abs( (h_enhancer_colapse$end + 1) - h_enhancer_colapse$start))
    
    # orthologous sequences for human annotation to get total alignment
    ann_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                      "chromosome_mouse", "start_mouse", "end_mouse")],
                                 bed = human_enhancer_sub[category == ann[i]][,c("chromosome", "start", "end")])
    n_aligned[i] <- sum(abs((ann_align$end_A + 1) - ann_align$start_A))
    
    # total of mouse annotation in orthologous sequences for conservation
    m_enhancer_conserved <- bed.intersect(bed1 = ann_align[,c("chromosome_B", "start_B", "end_B")],
                                     bed2 = mouse_enhancer_sub[category == ann[i]][,c("chromosome", "start", "end")])
    n_conserved[i] <- sum(abs( (m_enhancer_conserved$end + 1) - m_enhancer_conserved$start))
    
    print(ann[i])
  }
  
  out_list[[j]] <- data.table( tissue = names(human_enhancer_specific)[j],
                               chromosome = chr,
                               annotation = ann,
                               n_total = n_total,
                               n_aligned = n_aligned,
                               n_conserved = n_conserved)
  
  print(j)
}

out_dt_specific <- do.call("rbind", out_list)
out_dt_specific$p_aligned <- out_dt_specific$n_aligned / out_dt_specific$n_total
out_dt_specific$p_conserved <- out_dt_specific$n_conserved / out_dt_specific$n_total

### EXPORT

out_all <- rbind(out_dt_all, out_dt_2)
fwrite(out_all, out.all.file, append = T)
fwrite(out_dt_specific, out.specific.file, append = T)


#####
