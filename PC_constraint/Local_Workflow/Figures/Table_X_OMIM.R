### FIGURE 3 -- 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plotROC)
library(pROC)
library("cowplot")

### Test each gene set agianst orthologous genes that are not in any gene list. 


### IMPORT 

omim <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/OMIM_Petrovski2013_Dataset_S1.csv")
newC <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
RVIS_esp <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/RVIS_Petrovski2013_ESP.csv")
Lek <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")


### FORMAT

### # convert "" to NA
omim[omim==""] <- NA

# colnames
colnames(RVIS_esp) <- c("H_external_gene_name", "RVIS_ESP", "RVIS_ESP_percentile")
Lek <- Lek[,c("gene", "mis_z", "lof_z", "pLI")]
colnames(Lek) <- c("H_external_gene_name", "mis_z_Lek", "lof_z_Lek", "pLI_Lek")

df <- merge(Lek, RVIS_esp, all = T)
df <- merge(df, newC, all = T)
df <- df[!is.na(df$H_external_gene_name),]
df <- df[df$H_external_gene_name != "",]
df <- df[complete.cases(df),]

# Set gene lists
omim.all <- na.omit(unique(omim$`OMIM disease genes`))
non.OMIM <- df$H_external_gene_name
non.OMIM <- non.OMIM[!non.OMIM %in% omim.all]

recesive <- na.omit(unique(omim$`OMIM Recessive`))
df$OMIM.recessive <- NA
df$OMIM.recessive[df$H_external_gene_name %in% recesive] <- 1
df$OMIM.recessive[df$H_external_gene_name %in% non.OMIM] <- 0

haploinsufficient <- na.omit(unique(omim$`OMIM Haploinsufficiency`))
df$OMIM.hap <- NA
df$OMIM.hap[df$H_external_gene_name %in% haploinsufficient] <- 1
df$OMIM.hap[df$H_external_gene_name %in% non.OMIM] <- 0

dom.neg <- na.omit(unique(omim$`OMIM Dominant Negative`))
df$OMIM.domneg <- NA
df$OMIM.domneg[df$H_external_gene_name %in% dom.neg] <- 1
df$OMIM.domneg[df$H_external_gene_name %in% non.OMIM] <- 0

denovo <- na.omit(unique(omim$`OMIM de novo`))
df$OMIM.denovo <- NA
df$OMIM.denovo[df$H_external_gene_name %in% denovo] <- 1
df$OMIM.denovo[df$H_external_gene_name %in% non.OMIM] <- 0

omim.all <- na.omit(unique(omim$`OMIM disease genes`))
df$OMIM.all <- NA
df$OMIM.all[df$H_external_gene_name %in% omim.all] <- 1
df$OMIM.all[df$H_external_gene_name %in% non.OMIM] <- 0

# percentiles
percentile <- ecdf(df$RVIS_ESP)
df$RVIS_ESP_percentile <- percentile(df$RVIS_ESP)
percentile <- ecdf(df$mis_z_Lek)
df$mis_z_percentile_Lek <- percentile(df$mis_z_Lek)
percentile <- ecdf(df$pLI_Lek)
df$pLI_percentile_Lek <- percentile(df$pLI_Lek)
percentile <- ecdf(df$M_fun_Z_0.0001)
df$M_fun_Z_percentile_0.0001 <- percentile(df$M_fun_Z_0.0001)
percentile <- ecdf(df$H_fun_Z_0.001)
df$H_fun_Z_percentile_0.001 <- percentile(df$H_fun_Z_0.001)

### Table 1 -- GLM and ROC

kapow <- function(OMIM, VARS, df){
  
  x <- c(OMIM, VARS)
  data <- df[, x, with = F]
  
  col10 <- names(data)[-1]
  
  # glm.list <- vector("list", length(col10))
  Estimate.vec <- rep(NA, length(col10))
  Std_Err.vec <- rep(NA, length(col10))
  P.vec <- rep(NA, length(col10))
  ROC.vec <- rep(NA, length(col10))
  N.vec <- rep(NA, length(col10))
  
  for(i in seq_along(col10)){
    
    ind <- !is.na(data[,(i+1), with = F])[,1]
    glm.data <- data[ind,]
    mod <- glm(reformulate(col10[i], names(glm.data)[1]), data = glm.data, family=binomial)
    # summary(mod)
    
    N.vec[i] <- table(glm.data[,1])[2]
    Estimate.vec[i] <- summary(mod)$coefficients[2,1]
    Estimate.vec[i] <- round(Estimate.vec[i], digits = 3)
    Std_Err.vec[i] <- summary(mod)$coefficients[2,2]
    Std_Err.vec[i] <- round(Std_Err.vec[i], digits = 3)
    P.vec[i] <- summary(mod)$coefficients[2,4]
    P.vec[i] <- signif(P.vec[i], digits = 3)
  
    pred <- predict(mod, glm.data)
    obs <- as.numeric(unlist(glm.data[,1]))
    ROC.vec[i] <- auc(roc(obs, pred))
    ROC.vec[i] <- round(ROC.vec[i], digits = 3)
  }
  
  names(N.vec) <- col10
  names(Estimate.vec) <- col10
  names(Std_Err.vec) <- col10
  names(P.vec) <- col10
  names(ROC.vec) <- col10
  
  # Combine to dataframe
  out <- matrix(N.vec, nrow = 1)
  out <- rbind(out, Estimate.vec)
  out <- rbind(out, Std_Err.vec)
  out <- rbind(out, P.vec)
  out <- rbind(out, ROC.vec)
  rownames(out) <- paste(names(glm.data)[1], c("N", "Estimate", "Std Error", "P", "ROC"), sep = " ")
  
  return(out)
}


OMIM = c(
  "OMIM.all",
  "OMIM.denovo",
  "OMIM.domneg",
  "OMIM.hap",
  "OMIM.recessive"
)
VARS <- c(
  "mis_z_percentile_Lek",
  "pLI_Lek",
  "RVIS_ESP_percentile",
  "H_fun_Z_percentile_0.001",
  "M_fun_Z_percentile_0.0001"
)
out.list <- list()
for (i in 1:length(OMIM)){
  out.list[[i]] <- kapow(OMIM[i], VARS, df)
}
table1 <- do.call("rbind", out.list)
# table1 <- t(table1)
table1 <- as.data.frame(table1)
table1 <- cbind(rownames(table1), data.frame(table1, row.names=NULL))
# colnames(table1)[1] <- "Constraint metric"
fwrite(table1, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Tables/Table_OMIM_ROC.csv")

