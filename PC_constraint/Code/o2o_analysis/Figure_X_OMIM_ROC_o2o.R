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

# get o2o orthologues
newC <- newC[newC$orthology_type == "ortholog_one2one"]

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
    Std_Err.vec[i] <- summary(mod)$coefficients[2,2]
    P.vec[i] <- summary(mod)$coefficients[2,4]
    
    pred <- predict(mod, glm.data)
    obs <- as.numeric(unlist(glm.data[,1]))
    ROC.vec[i] <- auc(roc(obs, pred))
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

# remove NAs
table.df <- df
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
  out.list[[i]] <- kapow(OMIM[i], VARS, table.df)
}
table1 <- do.call("rbind", out.list)
# table1 <- t(table1)
# table1 <- as.data.frame(table1)
# table1 <- cbind(rownames(table1), data.frame(table1, row.names=NULL))
# colnames(table1)[1] <- "Constraint metric"
# fwrite(table1, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Tables/Table_OMIM_ROC.csv")

### Haploinsufficient ROC

OMIM <- "OMIM.hap"
kapow2 <- function(OMIM, VARS, df){
  
  x <- c(OMIM, VARS)
  glm.data <- df[, x, with = F]
  
  col10 <- names(glm.data)[-1]
  
  # glm.list <- vector("list", length(col10))
  Pred.list <- list()
  
  for(i in seq_along(col10)){
    mod <- glm(reformulate(col10[i], names(glm.data)[1]), data = glm.data, family=binomial)
    pred <- predict(mod, glm.data)
    obs <- as.numeric(unlist(glm.data[,1]))
    Pred.list[[i]] <- data.frame(H_external_gene_name = df$H_external_gene_name,
                                 obs = obs,
                                 pred = pred)
    colnames(Pred.list[[i]]) <- c("H_external_gene_name", "obs", paste(col10[i], colnames(Pred.list[[i]])[3], sep = "_"))
  }
  
  mod.all <- glm(reformulate(col10, OMIM), data = glm.data, family=binomial)
  summary(mod.all)
  pred <- predict(mod.all, glm.data)
  obs <- as.numeric(unlist(glm.data[,1]))
  auc(roc(obs, pred))
  Pred.list[[i+1]] <- data.frame(H_external_gene_name = df$H_external_gene_name,
                                 obs = obs,
                                 Combined_pred = pred)
  
  # names(Pred.list) <- c(col10, "Combined")
  
  return(Pred.list)
}

OMIM <- "OMIM.hap"
VARS <- c(
  "mis_z_percentile_Lek",
  "pLI_Lek",
  "RVIS_ESP_percentile",
  "H_fun_Z_percentile_0.001",
  "M_fun_Z_percentile_0.0001"
)
ROC_df <- kapow2(OMIM, VARS, table.df)
ROC_df <- do.call("cbind", ROC_df)
ROC_df <- ROC_df[, !duplicated(colnames(ROC_df))]


ROC_df1 <- ROC_df[,1:(ncol(ROC_df)-1)]

n.mizZ <- table1["OMIM.hap N", "mis_z_percentile_Lek"]
n.pLI <- table1["OMIM.hap N", "pLI_Lek"]
n.RVIS <- table1["OMIM.hap N", "RVIS_ESP_percentile"]
n.MfunZ <- table1["OMIM.hap N", "H_fun_Z_percentile_0.001"]
n.HfunZ <- table1["OMIM.hap N", "M_fun_Z_percentile_0.0001"]
auc.mizZ <- round(table1["OMIM.hap ROC", "mis_z_percentile_Lek"], 2)
auc.pLI <- round(table1["OMIM.hap ROC", "pLI_Lek"], 2)
auc.RVIS <- round(table1["OMIM.hap ROC", "RVIS_ESP_percentile"], 2)
auc.MfunZ <- round(table1["OMIM.hap ROC", "H_fun_Z_percentile_0.001"], 2)
auc.HfunZ <- round(table1["OMIM.hap ROC", "M_fun_Z_percentile_0.0001"], 2)
p.mizZ <- formatC(table1["OMIM.hap P", "mis_z_percentile_Lek"], format = "e", digits = 2)
p.pLI <- formatC(table1["OMIM.hap P", "pLI_Lek"], format = "e", digits = 2)
p.RVIS <- formatC(table1["OMIM.hap P", "RVIS_ESP_percentile"], format = "e", digits = 2)
p.MfunZ <- formatC(table1["OMIM.hap P", "H_fun_Z_percentile_0.001"], format = "e", digits = 2)
p.HfunZ <- formatC(table1["OMIM.hap P", "M_fun_Z_percentile_0.0001"], format = "e", digits = 2)

colnames(ROC_df1) <- c("H_external_gene_name", "obs", 
                       paste0("missense Z-score\nAUC=", auc.mizZ, " p=", p.mizZ),
                       paste0("pLI\nAUC=", auc.pLI, " p=", p.pLI),                  
                       paste0("RVIS\nAUC=", auc.RVIS, " p=", p.RVIS),    
                       paste0("Human funZ\nAUC=", auc.HfunZ, " p=", p.HfunZ),
                       paste0("Mouse funZ\nAUC=", auc.MfunZ, " p=", p.MfunZ))
ROC_df1 <- melt(ROC_df1, id.vars = c("H_external_gene_name", "obs"))
ROC_df1$variable <- as.factor(ROC_df1$variable)
unique(ROC_df1$variable)
ROC_df1$variable <- factor(ROC_df1$variable, levels = c(
  "pLI\nAUC=0.76 p=1.46e-25", 
  "Mouse funZ\nAUC=0.75 p=1.49e-22",
  "Human funZ\nAUC=0.75 p=1.02e-21",
  "missense Z-score\nAUC=0.74 p=8.10e-21", 
  "RVIS\nAUC=0.73 p=4.84e-20"
))
plotB <- ggplot(ROC_df1, aes(d = obs, m = value, color = variable)) + 
  geom_roc(n.cuts = 0) + 
  style_roc() +
  theme(legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.75, 0.25),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
plotB
# ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3B_o2o.jpg", plot = plotB, width = 6.5, height = 6)

### Multivartiant model 

multi_mod <- glm(OMIM.hap ~ pLI_Lek + 
                   M_fun_Z_percentile_0.0001 + 
                   H_fun_Z_percentile_0.001 +
                   mis_z_percentile_Lek +
                   RVIS_ESP_percentile, data = table.df, family = binomial)

pred <- predict(multi_mod, table.df)
obs <- as.numeric(table.df$OMIM.hap)
auc(roc(obs, pred))

p.mizZ <- formatC(summary(multi_mod)$coefficients["mis_z_percentile_Lek",4], format = "e", digits = 2)
p.pLI <- formatC(summary(multi_mod)$coefficients["pLI_Lek",4], format = "e", digits = 2)
p.RVIS <- formatC(summary(multi_mod)$coefficients["RVIS_ESP_percentile",4], format = "e", digits = 2)
p.MfunZ <- formatC(summary(multi_mod)$coefficients["M_fun_Z_percentile_0.0001",4], format = "e", digits = 2)
p.HfunZ <- formatC(summary(multi_mod)$coefficients["H_fun_Z_percentile_0.001",4], format = "e", digits = 2)

ROC_df2 <- ROC_df[,c("H_external_gene_name", "obs", "Combined_pred")]
colnames(ROC_df2) <- c("H_external_gene_name", "obs", 
                       paste0("Multivariate GLM AUC=0.80\n", 
                              "    pLI  p=", p.pLI, "\n",
                              "    Mouse funZ  p=", p.MfunZ, "\n",
                              "    RVIS  p=", p.RVIS, "\n",
                              "    missense Z  p=", p.mizZ, "\n",
                              "    Human funZ  p=", p.HfunZ))
ROC_df2 <- melt(ROC_df2, id.vars = c("H_external_gene_name", "obs"))
ROC_df2$variable <- as.factor(ROC_df2$variable)

plotC <- ggplot(ROC_df2, aes(d = obs, m = value, color = variable)) + 
  geom_roc(n.cuts = 0) + 
  style_roc() +
  scale_color_manual(values=c("black")) +
  theme(legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.65, 0.25),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
plotC
# ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3C_o2o.jpg", plot = plotC, width = 6.5, height = 6)


Fig3 <- plot_grid(plotB, plotC ,ncol = 2, nrow = 1, 
                  labels = c("A", "B"), label_size = 20)
Fig3
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3_o2o.jpg",
Fig3,
ncol = 1, nrow = 1,
base_height = 6, base_width = 13)
