rm(list = ls())
graphics.off()

library(data.table)

# Data: 
#   annotation
#   snvs
#   coverage_mask
#   soft mask
#   N mask
#   CGs
#   methylated CGs
# 
# covariates:
#   n SNVs
#   length
#   number of methylated CGs
#   mb SNV rate !?!?

### FUNCTIONS

# source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

bed_match <- function(bed, vec){
  
  require(data.table)
  
  bed_dt <- setDT(bed) # set as bed as data.table
  colnames(bed_dt) <- c("start", "end")
  bed_dt[, ind := .I] # add uniqe index to data.table
  setkey(bed_dt) # sets keys // order data by all columns
  
  vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
  vec_dt[, vec2 := vec] # dublicate column
  
  # Fast overlap join:
  ans1 = foverlaps(vec_dt, bed_dt, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
                   type = "within", nomatch = 0L)
  counts <- ans1[, .N, keyby = ind] # count by ind
  # merge to inital data
  bed_dt[, n_match := counts[bed_dt, on = .(ind), x.N]]
  setorder(bed_dt, ind) # reorder by ind to get inital order
  bed_dt[, ind := NULL] # deletes ind colum
  bed_dt[is.na(n_match), n_match := 0L] # NAs is 0 count
  
  return(bed_dt$n_match) 
}

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT 

ann <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz")
cg <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Methylation/Formatted/m.mus_GRC38_Ensembl_v101_CG_meth_NPC_ES.bed.gz")
sm <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_sm.bed.gz")
Nm <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz")

meth <- cg[percent_methylated > 60][,c("chromosome", "start", "end")]
cg <- cg[,c("chromosome", "start", "end")]

out_list <- list()
for (i in 1:19){
  
  snvs_chr <- fread(paste0("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/MGP/Variants/Formatted/MGP_v5_all_snps_PASS_chr", i, ".vcf.gz"))
  cm_chr <- fread(paste0("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/MGP/Coverage/Formatted/MGP_all_strain_DP_less10X_more0.1_chr", i, ".bed.gz"))
  ann_chr <- ann[chromosome == i]
  meth_chr <- meth[chromosome == i]
  cg_chr <- cg[chromosome == i]
  sm_chr <- sm[V1 == i]
  Nm_chr <- Nm[V1 == i]
  
  snvs_chr <- snvs_chr[, 1:2]
  snvs_chr$end <- snvs_chr$POS
  colnames(snvs_chr) <- c("chromosome", "start", "end")
  
  ann_chr$n_snv <- bed_match(bed = ann_chr[,c("start", "end")], 
                             vec = snvs_chr$start)
  
  ann_chr$n_Nmask <- bed_match(bed = ann_chr[,c("start", "end")], 
                               vec = unlist(seq2(from = Nm_chr$V2, to = Nm_chr$V3)))
  
  ann_chr$n_DPmask <- bed_match(bed = ann_chr[,c("start", "end")], 
                                vec = unlist(seq2(from = cm_chr$start, to = cm_chr$end)))
  
  ann_chr$n_soft_mask <- bed_match(bed = ann_chr[,c("start", "end")], 
                                   vec = unlist(seq2(from = sm_chr$V2, to = sm_chr$V3)))
  
  ann_chr$n_CG <- bed_match(bed = ann_chr[,c("start", "end")], 
                                 vec = cg_chr$start)
  
  ann_chr$n_meth_CG <- bed_match(bed = ann_chr[,c("start", "end")], 
                                 vec = meth_chr$start)
  
  ann_chr$length <- (ann_chr$end + 1) - ann_chr$start
  
  out_list[[i]] <- ann_chr
  print(i)
}

output <- do.call("rbind", out_list)

fwrite(output, 
       "~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_nonCDS_constraint/Sand/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell_covariates.csv.gz",
       compress = "gzip")

#####

# sum((cm$end + 1) - cm$start) 
# sum(cm$end[nrow(cm)] - cm$start[1])








