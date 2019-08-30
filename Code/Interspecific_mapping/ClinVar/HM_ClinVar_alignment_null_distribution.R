rm(list = ls())
graphics.off()

library(data.table)


### IMPORT

dt_list <- list()
for (i in 1:22){
  dt_list[[i]] <- fread(paste0("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_mapping_null_sampling_Hchr", i, ".csv"))
}
dt <- do.call("rbind", dt_list)
cv <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_mouse_mapping.csv")

### FORMAT 

colnames(dt) <- c("H_annotation", "N_SNVs", "N_align", "N_conserved", "N_conserved_K", "ind", "H_chromosome")
colnames(cv) <- c("H_CHR", "H_POS", "H_REF", "H_Kmer", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")

annotation <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

## sum null by annotation and int
ann_list <- list()
for (ann in 1:length(annotation)){
  sub <- subset(dt, dt$H_annotation %like% annotation[ann])
  ind_list <- list()
  for (i in 1:length(unique(sub$ind))){
    N_SNV <- sum(sub$N_SNVs[sub$ind == i])
    N_align <- sum(sub$N_align[sub$ind == i])
    N_conserved <- sum(sub$N_conserved[sub$ind == i])
    N_conserved_K <- sum(sub$N_conserved_K[sub$ind == i])
    ind_list[[i]] <- c(N_SNV, N_align, N_conserved, N_conserved_K)
    print(i)
  }
 tmp_ann <- as.data.table(do.call("rbind", ind_list))
 colnames(tmp_ann) <- c("N_SNV", "N_align", "N_conserved", "N_conserved_K")
 tmp_ann$H_annotation <- annotation[ann]
 tmp_ann$N_align_frac <- tmp_ann$N_align/tmp_ann$N_SNV # calculate fractions
 tmp_ann$N_conserved_frac <- tmp_ann$N_conserved/tmp_ann$N_SNV # calculate fractions
 tmp_ann$N_conserved_K_frac <- tmp_ann$N_conserved_K/tmp_ann$N_SNV # calculate fractions
 ann_list[[ann]] <- tmp_ann
}

## calculate null mean and sd by annotaion
ann_param <- list()
for (i in 1:length(ann_list)){
  N_SNV <- mean(ann_list[[i]]$N_SNV)  
  F_align_m <- mean(ann_list[[i]]$N_align_frac) 
  F_align_sd <- sd(ann_list[[i]]$N_align_frac) 
  F_conserved_m <- mean(ann_list[[i]]$N_conserved_frac) 
  F_conserved_sd <- sd(ann_list[[i]]$N_conserved_frac) 
  F_conserved_K_m <- mean(ann_list[[i]]$N_conserved_K_frac) 
  F_conserved_K_sd <- sd(ann_list[[i]]$N_conserved_K_frac) 
  ann_param[[i]] <- c(N_SNV, F_align_m, F_align_sd, F_conserved_m, F_conserved_sd, F_conserved_K_m, F_conserved_K_sd, annotation[i])
}
ann_param <- as.data.table(do.call("rbind", ann_param))
colnames(ann_param) <- c("N_SNV", "F_align_m", "F_align_sd", "F_conserved_m", "F_conserved_sd", "F_conserved_K_m", "F_conserved_K_sd", "H_annotation")

## sum ClinVar SNVs by annotation
cv_list <- list()
for (i in 1:length(annotation)){
  sub <- subset(cv, cv$H_ANN %like% annotation[i])
  N_SNV_obs <- nrow(sub)
  F_align_obs <- nrow(sub[!is.na(sub$M_POS),])/N_SNV_obs
  F_conserved_obs <- nrow(sub[sub$H_REF == sub$M_REF,])/N_SNV_obs
  F_conserved_K_obs <- nrow(sub[sub$H_Kmer == sub$M_Kmer,])/N_SNV_obs
  cv_list[[i]] <-  c(  annotation[i], N_SNV_obs, F_align_obs, F_conserved_obs, F_conserved_K_obs)
  names(cv_list[[i]]) <- c("H_annotation", "N_SNV_obs", "F_align_obs", "F_conserved_obs", "F_conserved_K_obs")
}
cv_obs <- as.data.table(do.call("rbind", cv_list))


## calculate Z scores 
out <- as.data.frame(cv_obs[ann_param, on = "H_annotation"])
cols_to_num <- c("N_SNV_obs", "F_align_obs", "F_conserved_obs", "F_conserved_K_obs", "N_SNV", "F_align_m", "F_align_sd",      
                 "F_conserved_m", "F_conserved_sd", "F_conserved_K_m", "F_conserved_K_sd")
out[cols_to_num] <- sapply(out[cols_to_num],as.numeric)
out$F_align_z <- (out$F_align_obs - out$F_align_m) / out$F_align_sd
out$F_conserved_z <- (out$F_conserved_obs - out$F_conserved_m) / out$F_conserved_sd
out$F_conserved_K_z <- (out$F_conserved_K_obs - out$F_conserved_K_m) / out$F_conserved_K_sd
# out <- melt(out, id.vars = "H_annotation")


### remove CDS

dt <- subset(dt, !dt$H_annotation %like% "A")
cv <- subset(cv, !cv$H_ANN %like% "A")
annotation <- c("B", "C", "D", "E", "F", "G", "H", "I", "J")

## sum null by annotation and int
ann_list <- list()
for (ann in 1:length(annotation)){
  sub <- subset(dt, dt$H_annotation %like% annotation[ann])
  ind_list <- list()
  for (i in 1:length(unique(sub$ind))){
    N_SNV <- sum(sub$N_SNVs[sub$ind == i])
    N_align <- sum(sub$N_align[sub$ind == i])
    N_conserved <- sum(sub$N_conserved[sub$ind == i])
    N_conserved_K <- sum(sub$N_conserved_K[sub$ind == i])
    ind_list[[i]] <- c(N_SNV, N_align, N_conserved, N_conserved_K)
    print(i)
  }
  tmp_ann <- as.data.table(do.call("rbind", ind_list))
  colnames(tmp_ann) <- c("N_SNV", "N_align", "N_conserved", "N_conserved_K")
  tmp_ann$H_annotation <- annotation[ann]
  tmp_ann$N_align_frac <- tmp_ann$N_align/tmp_ann$N_SNV # calculate fractions
  tmp_ann$N_conserved_frac <- tmp_ann$N_conserved/tmp_ann$N_SNV # calculate fractions
  tmp_ann$N_conserved_K_frac <- tmp_ann$N_conserved_K/tmp_ann$N_SNV # calculate fractions
  ann_list[[ann]] <- tmp_ann
}

## calculate null mean and sd by annotaion
ann_param <- list()
for (i in 1:length(ann_list)){
  N_SNV <- mean(ann_list[[i]]$N_SNV)  
  F_align_m <- mean(ann_list[[i]]$N_align_frac) 
  F_align_sd <- sd(ann_list[[i]]$N_align_frac) 
  F_conserved_m <- mean(ann_list[[i]]$N_conserved_frac) 
  F_conserved_sd <- sd(ann_list[[i]]$N_conserved_frac) 
  F_conserved_K_m <- mean(ann_list[[i]]$N_conserved_K_frac) 
  F_conserved_K_sd <- sd(ann_list[[i]]$N_conserved_K_frac) 
  ann_param[[i]] <- c(N_SNV, F_align_m, F_align_sd, F_conserved_m, F_conserved_sd, F_conserved_K_m, F_conserved_K_sd, annotation[i])
}
ann_param <- as.data.table(do.call("rbind", ann_param))
colnames(ann_param) <- c("N_SNV", "F_align_m", "F_align_sd", "F_conserved_m", "F_conserved_sd", "F_conserved_K_m", "F_conserved_K_sd", "H_annotation")

## sum ClinVar SNVs by annotation
cv_list <- list()
for (i in 1:length(annotation)){
  sub <- subset(cv, cv$H_ANN %like% annotation[i])
  N_SNV_obs <- nrow(sub)
  F_align_obs <- nrow(sub[!is.na(sub$M_POS),])/N_SNV_obs
  F_conserved_obs <- nrow(sub[sub$H_REF == sub$M_REF,])/N_SNV_obs
  F_conserved_K_obs <- nrow(sub[sub$H_Kmer == sub$M_Kmer,])/N_SNV_obs
  cv_list[[i]] <-  c(  annotation[i], N_SNV_obs, F_align_obs, F_conserved_obs, F_conserved_K_obs)
  names(cv_list[[i]]) <- c("H_annotation", "N_SNV_obs", "F_align_obs", "F_conserved_obs", "F_conserved_K_obs")
}
cv_obs <- as.data.table(do.call("rbind", cv_list))


## calculate Z scores 
out2 <- as.data.frame(cv_obs[ann_param, on = "H_annotation"])
cols_to_num <- c("N_SNV_obs", "F_align_obs", "F_conserved_obs", "F_conserved_K_obs", "N_SNV", "F_align_m", "F_align_sd",      
                 "F_conserved_m", "F_conserved_sd", "F_conserved_K_m", "F_conserved_K_sd")
out2[cols_to_num] <- sapply(out2[cols_to_num],as.numeric)
out2$F_align_z <- (out2$F_align_obs - out2$F_align_m) / out2$F_align_sd
out2$F_conserved_z <- (out2$F_conserved_obs - out2$F_conserved_m) / out2$F_conserved_sd
out2$F_conserved_K_z <- (out2$F_conserved_K_obs - out2$F_conserved_K_m) / out2$F_conserved_K_sd
# out <- melt(out, id.vars = "H_annotation")

### OUTPUT

fwrite(out2, "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_z_scores.csv")



