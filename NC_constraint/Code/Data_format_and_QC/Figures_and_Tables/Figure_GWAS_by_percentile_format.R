rm(list=ls())
graphics.off()

library(data.table)

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied", call.=FALSE)
# } 

H.gwas.file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_catalog_formatted.txt"
HM.gwas.file <- "~/Dropbox/PhD/Data/GWAS/formatted/GWAS_human_to_mouse_alignment.txt"
H.con.file <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_human_950_50_MAF001.csv"
M.con.file <- "~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_mouse_950_50.csv"
M.ann.file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Mouse_GRCh38_GENCODE_RegBuild_annotation_ranked.csv"
out.file <- "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/GWAS_SNVs_by_percentile.csv"


### FUNCTIONS
# seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


# IMPORT
H_CV <- fread(H.gwas.file)
HM_CV <- fread(HM.gwas.file)
H_con <- fread(H.con.file)
M_con <- fread(M.con.file)
M_ann <- fread(M.ann.file)

### FORMAT

# get alignment fraction for each percentile
M_con$length <- (M_con$POS_to + 1) - M_con$POS_from
M_frac <- merge(
  aggregate(n_aligned_central ~ Constraint_percentile_CpG, data=M_con, FUN=sum),
  aggregate(length ~ Constraint_percentile_CpG, data=M_con, FUN=sum)
)
M_frac$fraction <- M_frac$n_aligned_central / M_frac$length

H_con$length <- (H_con$POS_to + 1) - H_con$POS_from
H_frac <- merge(
  aggregate(n_aligned_central ~ Constraint_percentile_CpG, data=H_con, FUN=sum),
  aggregate(length ~ Constraint_percentile_CpG, data=H_con, FUN=sum)
)
H_frac$fraction <- H_frac$n_aligned_central / H_frac$length


H_CV <- H_CV[,c("CHR", "POS")]
# H_CV <- H_CV[H_CV$CHR %in% 1:22,]
H_CV <- H_CV[!duplicated(H_CV),]
H_CV$CHR <- as.integer(H_CV$CHR)
H_CV$POS <- as.integer(H_CV$POS)

HM_CV <- subset(HM_CV, HM_CV$M_CHR %in% 1:19)
HM_CV$M_CHR <- as.integer(HM_CV$M_CHR)
HM_CV_H <- HM_CV[,c("CHR","POS","REF")]
HM_CV_H <- setDT(HM_CV_H)[H_con, Constraint_percentile := Constraint_percentile_CpG, 
                          on = .(CHR, POS >= POS_from, POS <= POS_to)]
colnames(HM_CV_H) <- paste0("H_", colnames(HM_CV_H))
HM_CV_M <- HM_CV[,c("M_CHR","M_POS","M_REF")]
colnames(HM_CV_M) <- c("CHR","POS","REF")
HM_CV_M <- setDT(HM_CV_M)[M_con, Constraint_percentile := Constraint_percentile_CpG, 
                          on = .(CHR, POS >= POS_from, POS <= POS_to)]

colnames(M_ann) <- c("category","CHR","start","end") 
HM_CV_M <- setDT(HM_CV_M)[M_ann, CAT := category, 
                          on = .(CHR, POS >= start, POS <= end)]

colnames(HM_CV_M) <- paste0("M_", colnames(HM_CV_M))
HM_CV <- cbind(HM_CV_H, HM_CV_M)
HM_CV <- HM_CV[!duplicated(HM_CV),]

H_CV <- setDT(H_CV)[H_con, Constraint_percentile := Constraint_percentile_CpG, 
                    on = .(CHR, POS >= POS_from, POS <= POS_to)]
colnames(H_CV) <- paste0("H_", colnames(H_CV))

# summary(is.na(H_CV$H_Constraint_percentile))
# summary(is.na(HM_CV$H_Constraint_percentile))
# summary(is.na(HM_CV$M_Constraint_percentile))

dt <- HM_CV[H_CV, on = c("H_CHR", "H_POS", "H_Constraint_percentile"), allow.cartesian=TRUE]
dt <- dt[!duplicated(dt),]

summary(is.na(dt$H_Constraint_percentile))
summary(is.na(dt$M_Constraint_percentile))

# H_dt <- dt[,c("H_CHR", "H_POS", "H_Constraint_percentile", "H_CAT")]
# H_dt <- H_dt[complete.cases(H_dt),]
# H_dt$H_CAT[H_dt$H_CAT != "Exon - CDS"] <- "Non-coding"
# H_dt$H_CAT[H_dt$H_CAT == "Exon - CDS"] <- "Protein-coding"
# table(H_dt$H_CAT)
# 
# 
# M_dt <- dt[,c("M_CHR", "M_POS", "M_Constraint_percentile", "M_CAT")]
# M_dt <- M_dt[complete.cases(M_dt),]
# M_dt$M_CAT[M_dt$M_CAT != "Exon - CDS"] <- "Non-coding"
# M_dt$M_CAT[M_dt$M_CAT == "Exon - CDS"] <- "Protein-coding"
# table(M_dt$M_CAT)

percentile <- 1:100
H_n_GWAS <- rep(NA, 100)
M_n_GWAS <- rep(NA, 100)


for(i in percentile){
  tmp <- subset(dt, dt$H_Constraint_percentile == i)
  H_n_GWAS[i] <- nrow(tmp)
  tmp <- subset(dt, dt$M_Constraint_percentile == i)
  M_n_GWAS[i] <- nrow(tmp)
}

out <- data.frame(
                  percentile = c(percentile, percentile),
                  n_gwas = c(H_n_GWAS, M_n_GWAS),
                  fraction = c(rep(1, 100), M_frac$fraction),
                  species = c(rep("Human", 100), rep("Mouse", 100))
                  )
out$n_gwas_adj <- (out$n_gwas / out$fraction)



### EXPORT
fwrite(out, out.file)


# NC_dt <- subset(dt, dt$H_CAT == "Non-coding" & dt$M_CAT == "Non-coding")
# plot(NC_dt$H_Constraint_percentile, NC_dt$M_Constraint_percentile)
# cor.test(NC_dt$H_Constraint_percentile, NC_dt$M_Constraint_percentile)
# test <- aggregate(NC_dt$H_Constraint_percentile, by=list(Category=NC_dt$M_Constraint_percentile), FUN=median, na.rm = T)
# plot(test$x ~ test$Category)
# summary(lm(test$x ~ test$Category))
# 
# cor.test(dt$H_Constraint_percentile, dt$M_Constraint_percentile)
# test <- aggregate(dt$H_Constraint_percentile, by=list(Category=dt$M_Constraint_percentile), FUN=median, na.rm = T)
# plot(test$x ~ test$Category)
# summary(lm(test$x ~ test$Category))

#####

# # get vector of all annotation POS_ID 
# HM_CV <- subset(HM_CV, HM_CV$M_CHR != "X" & HM_CV$CHR != "X")
# H_POS_ID_all <- paste0(H_CV$CHR, "_", H_CV$POS)
# H_POS_ID_align <- paste0(HM_CV$CHR, "_", HM_CV$POS)
# M_POS_ID_align <- paste0(HM_CV$M_CHR, "_", HM_CV$M_POS)
# 
# 
# # get number of annotation POS in each percentile
# H_all <- rep(NA, 100)
# H_align <- rep(NA, 100)
# M_align <- rep(NA, 100)
# percentile <- c(1:100)
# for (i in percentile){
#   H.p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/human_POS_ID_percentile_", i, "_650_50.csv"), header = F)
#   M.p <- fread(paste0("/well/lindgren/George/Data/NC_constraint/Percentile_POS_ID/mouse_POS_ID_percentile_", i, "_650_50.csv"), header = F)
#   H_all[i] <- sum(H.p$V1 %chin% H_POS_ID_all)
#   H_align[i] <- sum(H.p$V1 %chin% H_POS_ID_align)
#   M_align[i] <- sum(M.p$V1 %chin% M_POS_ID_align)
#   print(i)
# }
# 
# df_out <- data.frame(Percentile = c(percentile, percentile, percentile), 
#                      N_gwas = c(H_all, H_align, M_align),
#                      Annotation = c(rep("Human gwas all", 100), rep("Human gwas aligned", 100), rep("Mouse gwas aligned", 100)))
