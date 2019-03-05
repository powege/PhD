### FIGURE 3 -- 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plotROC)
library(pROC)
library("cowplot")

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

# percentiles
RVIS_esp$RVIS_ESP_percentile <- RVIS_esp$RVIS_ESP_percentile/100
percentile <- ecdf(Lek$mis_z_Lek)
Lek$mis_z_percentile_Lek <- percentile(Lek$mis_z_Lek)
percentile <- ecdf(Lek$pLI_Lek)
Lek$pLI_percentile_Lek <- percentile(Lek$pLI_Lek)

df <- merge(Lek, RVIS_esp, all = T)
df <- merge(df, newC, all = T)
df <- df[!is.na(df$H_external_gene_name),]
df <- df[df$H_external_gene_name != "",]

# Set gene lists
recesive <- na.omit(unique(omim$`OMIM Recessive`))
df$OMIM.recessive <- 0
df$OMIM.recessive[df$H_external_gene_name %in% recesive] <- 1

haploinsufficient <- na.omit(unique(omim$`OMIM Haploinsufficiency`))
df$OMIM.hap <- 0
df$OMIM.hap[df$H_external_gene_name %in% haploinsufficient] <- 1

dom.neg <- na.omit(unique(omim$`OMIM Dominant Negative`))
df$OMIM.domneg <- 0
df$OMIM.domneg[df$H_external_gene_name %in% dom.neg] <- 1

denovo <- na.omit(unique(omim$`OMIM de novo`))
df$OMIM.denovo <- 0
df$OMIM.denovo[df$H_external_gene_name %in% denovo] <- 1

omim.all <- na.omit(unique(omim$`OMIM disease genes`))
df$OMIM.all <- 0
df$OMIM.all[df$H_external_gene_name %in% omim.all] <- 1

# subset one to one orthologs
# df <- subset(df, df$orthology_type == "ortholog_one2one")

### N OMIM groups 
# n.recessive <- length(df$OMIM.recessive[df$OMIM.recessive == 1])
# n.hap <- length(df$OMIM.hap[df$OMIM.hap == 1])
# n.domneg <- length(df$OMIM.domneg[df$OMIM.domneg == 1])
# n.denovo <- length(df$OMIM.denovo[df$OMIM.denovo == 1])
# n.all <- length(df$OMIM.all[df$OMIM.all == 1])


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
table.df <- df[complete.cases(df), ]
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



### OMIM Boxplot
# x <-c("M_fun_Z_percentile_0.0001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all")
boxplot.OMIM <- function(x, df){
  
  box.df <- df[ ,x, with = F]
  colnames(box.df) <- c(x[1], 
                        paste0("Recessive"),
                        paste0("Haploinsufficient"),
                        paste0("Dominant\nnegative"),
                        paste0("Denovo"),
                        paste0("All"))
  box.df <- melt(box.df, id.vars = x[1])
  box.df$value[box.df$value == 1] <- "In OMIM set"
  box.df$value[box.df$value == 0] <- "Not in OMIM set"
  ind <- which(box.df[,1] > -2 & box.df[,1] < 2)
  box.df <- box.df[ind,]
  
  box.df <- as.data.frame(box.df)
  box.df$value <- as.factor(box.df$value)
  out.box <- ggplot(box.df, aes_string(x=box.df[,2], y=box.df[,1], fill=box.df[,3])) +
    geom_boxplot() +
    ylab("Mouse funZ percentile") +
    xlab("") +
    scale_fill_manual(values=c("gray70", "gray100", "gray100")) +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text=element_text(size=14),
          legend.position="top",
          plot.margin=unit(c(1,1,0,1),"cm")) 
  out.box
  # medians <- box.df[, .(median(value)), by = c("IMPC.Viability.2", "variable")][order(IMPC.Viability.2, variable)]
  return(out.box)
}

plotA <- boxplot.OMIM(c("M_fun_Z_percentile_0.0001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"), df)
plotA

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
ROC_df <- kapow2(OMIM, VARS, df)
ROC_df <- do.call("cbind", ROC_df)
ROC_df <- ROC_df[, !duplicated(colnames(ROC_df))]
colnames(ROC_df) <- c("H_external_gene_name", "obs", 
                      "missense Z-score\n(AUC = 0.747)", 
                      "pLI***\n(AUC = 0.782)",                  
                      "RVIS*\n(AUC = 0.735)",    
                      "Human funZ\n(AUC = 0.759)",
                      "Mouse funZ***\n(AUC = 0.762)",
                      "Multivariate\n(AUC = 0.817)")
ROC_df <- melt(ROC_df, id.vars = c("H_external_gene_name", "obs"))
ROC_df$variable <- as.factor(ROC_df$variable)
ROC_df$variable <- factor(ROC_df$variable, levels = c("Multivariate\n(AUC = 0.817)",
                                                      "pLI***\n(AUC = 0.782)", 
                                                      "Mouse funZ***\n(AUC = 0.762)",
                                                      "Human funZ\n(AUC = 0.759)",
                                                      "missense Z-score\n(AUC = 0.747)", 
                                                      "RVIS*\n(AUC = 0.735)"))
plotB <- ggplot(ROC_df, aes(d = obs, m = value, color = variable)) + 
  geom_roc(n.cuts = 0) + 
  style_roc() +
  theme(legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
plotB


Fig3 <- plot_grid(plotA, plotB ,ncol = 1, nrow = 2, 
                  labels = c("A", "B"), label_size = 20)
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3.jpg", 
          Fig3, 
          ncol = 1, nrow = 1, 
          base_height = 10, base_width = 8)
ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3A.jpg", plot = plotA, width = 8, height = 6)
ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3B.jpg", plot = plotB, width = 8, height = 6)
 
