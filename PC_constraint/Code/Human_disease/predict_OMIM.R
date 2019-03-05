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



#################

# rm(list = ls())
# graphics.off()
# 
# library(data.table)
# library(pROC)
# library(MASS)
# 
# 
# ### IMPORT 
# 
# omim <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/OMIM_Petrovski2013_Dataset_S1.csv")
# newC <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
# RVIS_esp <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/RVIS_Petrovski2013_ESP.csv")
# Lek <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")
# 
# ### FORMAT
# 
# ### # convert "" to NA
# omim[omim==""] <- NA
# 
# # colnames
# colnames(RVIS_esp) <- c("H_external_gene_name", "RVIS_ESP", "RVIS_ESP_percentile")
# Lek <- Lek[,c("gene", "mis_z", "lof_z", "pLI")]
# colnames(Lek) <- c("H_external_gene_name", "mis_z_Lek", "lof_z_Lek", "pLI_Lek")
# 
# 
# # percentiles
# RVIS_esp$RVIS_ESP_percentile <- RVIS_esp$RVIS_ESP_percentile/100
# percentile <- ecdf(Lek$mis_z_Lek)
# Lek$mis_z_percentile_Lek <- percentile(Lek$mis_z_Lek)
# percentile <- ecdf(Lek$lof_z_Lek)
# Lek$lof_z_percentile_Lek <- percentile(Lek$lof_z_Lek)
# percentile <- ecdf(Lek$pLI_Lek)
# Lek$pLI_percentile_Lek <- percentile(Lek$pLI_Lek)
# 
# df <- merge(Lek, RVIS_esp, all = T)
# df <- merge(df, newC, all = T)
# df <- df[!is.na(df$H_external_gene_name),]
# df <- df[df$H_external_gene_name != "",]
# 
# # Set gene lists
# recesive <- na.omit(unique(omim$`OMIM Recessive`))
# df$OMIM.recessive <- 0
# df$OMIM.recessive[df$H_external_gene_name %in% recesive] <- 1
# 
# haploinsufficient <- na.omit(unique(omim$`OMIM Haploinsufficiency`))
# df$OMIM.hap <- 0
# df$OMIM.hap[df$H_external_gene_name %in% haploinsufficient] <- 1
# 
# dom.neg <- na.omit(unique(omim$`OMIM Dominant Negative`))
# df$OMIM.domneg <- 0
# df$OMIM.domneg[df$H_external_gene_name %in% dom.neg] <- 1
# 
# denovo <- na.omit(unique(omim$`OMIM de novo`))
# df$OMIM.denovo <- 0
# df$OMIM.denovo[df$H_external_gene_name %in% denovo] <- 1
# 
# omim.all <- na.omit(unique(omim$`OMIM disease genes`))
# df$OMIM.all <- 0
# df$OMIM.all[df$H_external_gene_name %in% omim.all] <- 1
# 
# ### N OMIM groups 
# n.recessive <- length(df$OMIM.recessive[df$OMIM.recessive == 1])
# n.hap <- length(df$OMIM.hap[df$OMIM.hap == 1])
# n.domneg <- length(df$OMIM.domneg[df$OMIM.domneg == 1])
# n.denovo <- length(df$OMIM.denovo[df$OMIM.denovo == 1])
# n.all <- length(df$OMIM.all[df$OMIM.all == 1])
# 
# 
# # x <-c("H_funZ_percentile", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all")
# boxplot.OMIM <- function(x){
#   
#   box.df <- df[ ,x, with = F]
#   colnames(box.df) <- c(x[1], 
#                         paste0("Rec"),
#                         paste0("Hap"),
#                         paste0("DNeg"),
#                         paste0("Denovo"),
#                         paste0("All"))
#   box.df <- melt(box.df, id.vars = x[1])
#   box.df$value[box.df$value == 1] <- "In OMIM set"
#   box.df$value[box.df$value == 0] <- "Not in OMIM set"
#   ind <- which(box.df[,1] > -2 & box.df[,1] < 2)
#   box.df <- box.df[ind,]
#   
#   box.df <- as.data.frame(box.df)
#   box.df$value <- as.factor(box.df$value)
#   out.box <- ggplot(box.df, aes_string(x=box.df[,2], y=box.df[,1], fill=box.df[,3])) +
#     geom_boxplot() +
#     ylab("Percentile") +
#     xlab("") +
#     scale_fill_manual(values=c("gray70", "gray100", "gray100")) +
#     theme_classic() +
#     theme(legend.title=element_blank(),
#           axis.text = element_text(size = 14),
#           axis.title = element_text(size = 14),
#           legend.text=element_text(size=14),
#           legend.position="top") 
#   out.box
#   # medians <- box.df[, .(median(value)), by = c("IMPC.Viability.2", "variable")][order(IMPC.Viability.2, variable)]
#   return(out.box)
# }
# 
# # boxplot.OMIM(c("H_fun_Z_percentile_0.001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# # boxplot.OMIM(c("H_fun_Z_percentile_0.0005", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# # 
# # boxplot.OMIM(c("H_mis_Z_percentile_0.001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# # boxplot.OMIM(c("H_mis_Z_percentile_0.0005", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# # boxplot.OMIM(c("H_mis_Z_percentile_0.0001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# # 
# # boxplot.OMIM(c("H_RVIS_percentile_0.001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# # boxplot.OMIM(c("H_RVIS_percentile_0.0005", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"))
# 
# 
# ### GLM
# 
# # x <- c("OMIM.denovo", 
# #        "mis_z_percentile_Lek","lof_z_percentile_Lek",
# #        "pLI_Lek","RVIS_ESP_percentile",
# #        "H_fun_Z_0.001", "H_fun_Z_0.0005", "H_fun_Z_0.0001",
# #        "H_mis_Z_0.001", "H_mis_Z_0.0005", "H_mis_Z_0.0001",
# #        "H_RVIS_0.001", "H_RVIS_0.0005",
# #        "M_RVIS_0.0001", "M_mis_Z_0.0001", "M_fun_Z_0.0001")
# kapow <- function(x){
#   
#   glm.data <- df[, x, with = F]
#   col10 <- names(glm.data)[-1]
#   
#   glm.list <- vector("list", length(col10))
#   ROC.vec <- rep(NA, length(col10))
#   
#   for(i in seq_along(col10)){
#     glm.list[[i]] <- glm(reformulate(col10[i], names(glm.data)[1]), data = glm.data, family=binomial)
#     pred <- predict(glm.list[[i]], glm.data)
#     obs <- as.numeric(unlist(glm.data[,1]))
#     ROC.vec[i] <- auc(roc(obs, pred))
#   }
#   
#   names(ROC.vec) <- col10
#   out <- list(glm.list, ROC.vec)
#   return(out)
# }
# 
# OMIM.rec <- kapow(c("OMIM.recessive", 
#                     "mis_z_percentile_Lek","lof_z_percentile_Lek",
#                     "pLI_Lek","RVIS_ESP_percentile",
#                     "H_RVIS_percentile_0.001",
#                     "H_fun_Z_percentile_0.001", 
#                     "H_mis_Z_percentile_0.0005",
#                     "M_RVIS_percentile_0.0001", 
#                     "M_fun_Z_percentile_0.0001",
#                     "M_mis_Z_percentile_0.0001"))[[2]]
# OMIM.hap <- kapow(c("OMIM.hap", 
#                     "mis_z_percentile_Lek","lof_z_percentile_Lek",
#                     "pLI_Lek","RVIS_ESP_percentile",
#                     "H_RVIS_percentile_0.001",
#                     "H_fun_Z_percentile_0.001", 
#                     "H_mis_Z_percentile_0.0005",
#                     "M_RVIS_percentile_0.0001", 
#                     "M_fun_Z_percentile_0.0001",
#                     "M_mis_Z_percentile_0.0001"))[[2]]
# OMIM.domneg <- kapow(c("OMIM.domneg", 
#                        "mis_z_percentile_Lek","lof_z_percentile_Lek",
#                        "pLI_Lek","RVIS_ESP_percentile",
#                        "H_RVIS_percentile_0.001",
#                        "H_fun_Z_percentile_0.001", 
#                        "H_mis_Z_percentile_0.0005",
#                        "M_RVIS_percentile_0.0001", 
#                        "M_fun_Z_percentile_0.0001",
#                        "M_mis_Z_percentile_0.0001"))[[2]]
# OMIM.denovo <- kapow(c("OMIM.denovo", 
#                        "mis_z_percentile_Lek","lof_z_percentile_Lek",
#                        "pLI_Lek","RVIS_ESP_percentile",
#                        "H_RVIS_percentile_0.001",
#                        "H_fun_Z_percentile_0.001", 
#                        "H_mis_Z_percentile_0.0005",
#                        "M_RVIS_percentile_0.0001", 
#                        "M_fun_Z_percentile_0.0001",
#                        "M_mis_Z_percentile_0.0001"))[[2]]
# OMIM.all <- kapow(c("OMIM.all", 
#                     "mis_z_percentile_Lek","lof_z_percentile_Lek",
#                     "pLI_Lek","RVIS_ESP_percentile",
#                     "H_RVIS_percentile_0.001",
#                     "H_fun_Z_percentile_0.001", 
#                     "H_mis_Z_percentile_0.0005",
#                     "M_RVIS_percentile_0.0001", 
#                     "M_fun_Z_percentile_0.0001",
#                     "M_mis_Z_percentile_0.0001"))[[2]]
# 
# # Combine to dataframe
# test <- matrix(OMIM.all, nrow = 1)
# rownames(test) <- "OMIM.all"
# test <- rbind(test, OMIM.denovo)
# test <- rbind(test, OMIM.domneg)
# test <- rbind(test, OMIM.hap)
# test <- rbind(test, OMIM.rec)
# 
# 
# 
# #####
# OMIM = c("OMIM.all",
#                          "OMIM.denovo",
#                          "OMIM.domneg",
#                          "OMIM.hap",
#                          "OMIM.recessive")
# Vars = c(
#        # "mis_z_percentile_Lek",
#        # "lof_z_percentile_Lek",
#        "pLI_Lek",
#        # "RVIS_ESP_percentile",
#        "H_fun_Z_percentile_0.001",
#        # "H_RVIS_percentile_0.001",
#        # "M_RVIS_percentile_0.0001",
#        "M_fun_Z_percentile_0.0001"
#        )
# kapow2 <- function(OMIM, VARS){
#   
#   x <- c(OMIM, Vars)
#   glm.data <- df[, x, with = F]
#   col10 <- names(glm.data)[(length(OMIM) +1) : ncol(glm.data)]
#   
#   glm.list <- vector("list", length(OMIM))
#   coef.list <- vector("list", length(OMIM))
#   ROC.vec <- rep(NA, length(col10))
#   
#   for(i in seq_along(OMIM)){
#     
#     # glm.list[[i]] <- glm(reformulate(col10, OMIM[i]), data = glm.data, family=binomial)
#     glm.list[[i]] <- glm(reformulate(col10, OMIM[i]), data = glm.data, family=binomial)
#     coef.list[[i]] <- coef(summary(glm.list[[i]]))
#   
#     pred <- predict(glm.list[[i]], glm.data)
#     obs <- as.numeric(unlist(glm.data[, OMIM[i], with = F]))
#     ROC.vec[i] <- auc(roc(obs, pred))
#   }
#   
#   names(ROC.vec) <- OMIM
#   out <- list(glm.list, coef.list, ROC.vec)
#   return(out)
# }
# 
# test2 <- kapow2(OMIM, VARS)
# test2[[3]]
# test2[[2]]
# 
# #####
# 
# # OMIM.rec <- kapow(c("OMIM.recessive", 
# #                     "mis_z_percentile_Lek","lof_z_percentile_Lek",
# #                     "pLI_Lek","RVIS_ESP_percentile",
# #                     "H_fun_Z_0.001", "H_fun_Z_0.0005", "H_fun_Z_0.0001",
# #                     "H_mis_Z_0.001", "H_mis_Z_0.0005", "H_mis_Z_0.0001",
# #                     "H_RVIS_0.001", "H_RVIS_0.0005",
# #                     "M_RVIS_0.0001", "M_mis_Z_0.0001", "M_fun_Z_0.0001",
# #                     "HM_mis_Z_0.001", "HM_RVIS_0.001"))[[2]]
# # OMIM.hap <- kapow(c("OMIM.hap", 
# #                     "mis_z_percentile_Lek","lof_z_percentile_Lek",
# #                     "pLI_Lek","RVIS_ESP_percentile",
# #                     "H_fun_Z_0.001", "H_fun_Z_0.0005", "H_fun_Z_0.0001",
# #                     "H_mis_Z_0.001", "H_mis_Z_0.0005", "H_mis_Z_0.0001",
# #                     "H_RVIS_0.001", "H_RVIS_0.0005",
# #                     "M_RVIS_0.0001", "M_mis_Z_0.0001", "M_fun_Z_0.0001",
# #                     "HM_mis_Z_0.001", "HM_RVIS_0.001"))[[2]]
# # OMIM.domneg <- kapow(c("OMIM.domneg", 
# #                        "mis_z_percentile_Lek","lof_z_percentile_Lek",
# #                        "pLI_Lek","RVIS_ESP_percentile",
# #                        "H_fun_Z_0.001", "H_fun_Z_0.0005", "H_fun_Z_0.0001",
# #                        "H_mis_Z_0.001", "H_mis_Z_0.0005", "H_mis_Z_0.0001",
# #                        "H_RVIS_0.001", "H_RVIS_0.0005",
# #                        "M_RVIS_0.0001", "M_mis_Z_0.0001", "M_fun_Z_0.0001",
# #                        "HM_mis_Z_0.001", "HM_RVIS_0.001"))[[2]]
# # OMIM.denovo <- kapow(c("OMIM.denovo", 
# #                        "mis_z_percentile_Lek","lof_z_percentile_Lek",
# #                        "pLI_Lek","RVIS_ESP_percentile",
# #                        "H_fun_Z_0.001", "H_fun_Z_0.0005", "H_fun_Z_0.0001",
# #                        "H_mis_Z_0.001", "H_mis_Z_0.0005", "H_mis_Z_0.0001",
# #                        "H_RVIS_0.001", "H_RVIS_0.0005",
# #                        "M_RVIS_0.0001", "M_mis_Z_0.0001", "M_fun_Z_0.0001",
# #                        "HM_mis_Z_0.001", "HM_RVIS_0.001"))[[2]]
# # OMIM.all <- kapow(c("OMIM.all", 
# #                     "mis_z_percentile_Lek","lof_z_percentile_Lek",
# #                     "pLI_Lek","RVIS_ESP_percentile",
# #                     "H_fun_Z_0.001", "H_fun_Z_0.0005", "H_fun_Z_0.0001",
# #                     "H_mis_Z_0.001", "H_mis_Z_0.0005", "H_mis_Z_0.0001",
# #                     "H_RVIS_0.001", "H_RVIS_0.0005",
# #                     "M_RVIS_0.0001", "M_mis_Z_0.0001", "M_fun_Z_0.0001",
# #                     "HM_mis_Z_0.001", "HM_RVIS_0.001"))[[2]]
# 
