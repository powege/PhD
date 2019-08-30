rm(list = ls())
graphics.off()

library(data.table)
library(plotROC)
library(pROC)
library("cowplot")

# Table -- comparing ROC between constraint metrics on gene lists
# Fig A -- boxplots of gene lists (dif from mean not in list)
# Fig B -- ROC plot of gene lists (in gene list vs not in gene list)

### Script that:
  # Formats list of genes
  # Formats constraint metrics
  # Creates table
  # Creates figure


##############
### IMPORT RAW
##############

## Gene lists 
# hap <- fread("~/Dropbox/PhD/Data/OMIM/OMIM-Gene-Map-Retrieval_haploinsufficiency_210819.csv")
# denovo <- fread("~/Dropbox/PhD/Data/OMIM/OMIM-Gene-Map-Retrieval_de-novo_210819.csv")
# domneg <- fread("~/Dropbox/PhD/Data/OMIM/OMIM-Gene-Map-Retrieval_dominant-negative_210819.csv")
# rec <- fread("~/Dropbox/PhD/Data/OMIM/OMIM-Gene-Map-Retrieval_recessive_210819.csv")
drug <- fread("~/Dropbox/PhD/Data/Drug_Targets/Drug_targets_list_for_George_P.csv", fill = T, header = F)
omim <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/OMIM_Petrovski2013_Dataset_S1.csv")

### Constraint scores
funz <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
rvis <- fread("~/Dropbox/PhD/Data/Constraint/RVIS_Petrovski2013_ESP.csv")
gnomad <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt")

##############################
### FORMAT constraint metrics
#############################

gnomad <- gnomad[,c("gene", "mis_z", "pLI", "oe_lof_upper")]
colnames(gnomad) <- c("H_external_gene_name", "mis_z", "pLI", "oe_lof_upper")
rvis <- rvis[,c("HGNC gene", "Residual Variation Intolerance Score")]
colnames(rvis) <- c("H_external_gene_name", "RVIS")
funz <- funz[,c("H_external_gene_name", "M_external_gene_name", "M_fun_Z", "H_fun_Z", "orthology_type")]

# merge
metrics <- gnomad[rvis, on = "H_external_gene_name"]
metrics <- metrics[funz, on = "H_external_gene_name"]
metrics <- metrics[!duplicated(metrics),]
# metrics <- metrics[complete.cases(metrics),]

# calculate percentiles
percentile <- ecdf(metrics$mis_z[!duplicated(metrics$H_external_gene_name)])
metrics$mis_z_percentile <- percentile(metrics$mis_z)
percentile <- ecdf(metrics$pLI[!duplicated(metrics$H_external_gene_name)])
metrics$pLI_percentile <- percentile(metrics$pLI)
percentile <- ecdf(metrics$oe_lof_upper[!duplicated(metrics$H_external_gene_name)])
metrics$oe_lof_upper_percentile <- percentile(metrics$oe_lof_upper)
percentile <- ecdf(metrics$RVIS[!duplicated(metrics$H_external_gene_name)])
metrics$RVIS_percentile <- percentile(metrics$RVIS)
percentile <- ecdf(metrics$M_fun_Z[!duplicated(metrics$M_external_gene_name)])
metrics$M_fun_Z_percentile <- percentile(metrics$M_fun_Z)
percentile <- ecdf(metrics$H_fun_Z[!duplicated(metrics$H_external_gene_name)])
metrics$H_fun_Z_percentile <- percentile(metrics$H_fun_Z)

# one_to_one orthologues
# metrics <- subset(metrics, metrics$orthology_type == "ortholog_one2one")

######################
### FORMAT gene lists
#####################

## Identify genes in each gene list
rec_gene <- na.omit(unique(omim$`OMIM Recessive`))
hap_gene <- na.omit(unique(omim$`OMIM Haploinsufficiency`))
domneg_gene <- na.omit(unique(omim$`OMIM Dominant Negative`))
denovo_gene <- na.omit(unique(omim$`OMIM de novo`))

length(rec_gene)
length(hap_gene)
length(domneg_gene)
length(denovo_gene)

# # vec <- hap$`Gene/Locus`
gene_split <- function(vec){
  out <- unlist(strsplit(vec, split=","))
  out <- gsub(" ", "", out)
  out <- sort(unique(out))
  return(out)
}
# hap_gene <- gene_split(hap$`Gene/Locus`)
# denovo_gene <- gene_split(denovo$`Gene/Locus`)
# domneg_gene <- gene_split(domneg$`Gene/Locus`)
# rec_gene <- gene_split(rec$`Gene/Locus`)
drug_gene <- gene_split(drug$V1)

# List of gene lists
gene_list <- list(rec_gene, hap_gene, domneg_gene, denovo_gene)
gene_list_names <- c("OMIM recessive", "OMIM hoploinsufficiency", "OMIM dominant-negative", "OMIM de novo")

## Identify genes not in each gene list
gene_list_neg <- list()
for (i in 1:length(gene_list)){
  gene_list_neg[[i]] <- unique(metrics$H_external_gene_name[which(!metrics$H_external_gene_name %in% gene_list[[i]])])
}

gene_list_dt <- list()
for (i in 1:length(gene_list)){
  gene_list_dt[[i]] <- data.table(H_external_gene_name = c(gene_list[[i]], gene_list_neg[[i]]),
                       gene_list = rep(gene_list_names[i], length(c(gene_list[[i]], gene_list_neg[[i]]))),
                       sig = c(rep(1, length(gene_list[[i]])), rep(0, length(gene_list_neg[[i]]))))
  }
gene_list_dt <- do.call("rbind", gene_list_dt)
gene_list_dt <- subset(gene_list_dt, gene_list_dt$H_external_gene_name != "")


##############
### WORKING DT
##############

dt <- gene_list_dt[metrics, on = "H_external_gene_name"]
dt <- dt[!is.na(dt$gene_list),]
dt$gene_list <- as.factor(dt$gene_list)
dt_long <- dt[,c("H_external_gene_name", "gene_list", "sig", "mis_z_percentile",      
                 "pLI_percentile", "oe_lof_upper_percentile", "RVIS_percentile",      
                 "M_fun_Z_percentile", "H_fun_Z_percentile")] 
dt_long <- melt(dt_long, 
                id.vars = c("H_external_gene_name","gene_list","sig"), 
                measure.vars = c("mis_z_percentile", "pLI_percentile", 
                                 "oe_lof_upper_percentile", "RVIS_percentile", 
                                 "M_fun_Z_percentile", "H_fun_Z_percentile"))

#############
### TABLE ROC 
############

scores <- levels(dt_long$variable)
terms <- levels(dt_long$gene_list)
terms_ls <- list()

for (i in 1:length(terms)){
  
  sub_term <- subset(dt_long, dt_long$gene_list == terms[i])
  
  N_sig <- rep(NA, length(scores))
  N_nosig <- rep(NA, length(scores))
  mean_sig <- rep(NA, length(scores))
  mean_nosig <- rep(NA, length(scores))
  estimate <- rep(NA, length(scores))
  std_er <- rep(NA, length(scores))
  p_val <- rep(NA, length(scores))
  ROC <- rep(NA, length(scores))
  
  for (j in 1:length(scores)){
    sub_score <- subset(sub_term, sub_term$variable == scores[j])
    sub_score <- sub_score[!duplicated(sub_score)]
    
    N_sig[j] <- length(sub_score$H_external_gene_name[sub_score$sig == 1 & !is.na(sub_score$value)])
    N_nosig[j] <- length(sub_score$H_external_gene_name[sub_score$sig == 0 & !is.na(sub_score$value)])
    
    mean_sig[j] <- mean(sub_score$value[sub_score$sig == 1], na.rm = T)
    mean_nosig[j] <- mean(sub_score$value[sub_score$sig == 0], na.rm = T)
    
    mod <- glm(sig ~ value, data = sub_score, family=binomial)
    
    estimate[j] <- summary(mod)$coefficients[2,1]
    std_er[j] <- summary(mod)$coefficients[2,2]
    p_val[j] <- summary(mod)$coefficients[2,4]
    
    pred <- predict(mod, sub_score)
    obs <- as.numeric(unlist(sub_score$sig))
    ROC[j] <- auc(roc(obs, pred))
    
  }
  
  terms_ls[[i]] <- data.frame(
    gene_list = rep(terms[i], length(scores)),
    constraint_score = scores,
    n_gene_list = N_sig,
    n_not_in_gene_list = N_nosig,
    mean_gene_list = round(mean_sig, digits = 3),
    mean_not_in_gene_list = round(mean_nosig, digits = 3),
    estimate = round(estimate, digits = 3),
    standard_error = round(std_er, digits = 3),
    p_val = formatC(p.adjust(p_val, method = "bonferroni", n = length(p_val)), format = "e", digits = 2),
    ROC = round(ROC, digits = 2))

}

ROC_table <- do.call("rbind", terms_ls)

##################
### FIGURE boxplot
##################

dt_MfunZ <- subset(dt_long, dt_long$variable == "M_fun_Z_percentile")

# calculate mean funZ percentile for each null group
null_average <- subset(dt_MfunZ, dt_MfunZ$sig == 0)
null_average <- aggregate(value ~ gene_list, data=null_average, FUN=mean)
colnames(null_average) <- c("gene_list", "null_average")

# calculate the relative difference in funZ for each group relative to null average
dt_plot1 <- subset(dt_MfunZ, dt_MfunZ$sig == 1)
dt_plot1 <- merge(dt_plot1, null_average, by = "gene_list")
dt_plot1$dif <- dt_plot1$value - dt_plot1$null_average

# set order of top-level terms
tmp1 <- subset(dt_MfunZ, dt_MfunZ$sig == 1)
order1 <- aggregate(tmp1$value~tmp1$gene_list, FUN=median)
tmp2 <- subset(dt_MfunZ, dt_MfunZ$sig == 0)
order2 <- aggregate(tmp2$value~tmp2$gene_list, FUN=mean)
colnames(order1) <- c("gene_list", "M_fun_Z_percentile_1")
colnames(order2) <- c("gene_list", "M_fun_Z_percentile_2")
order <- merge(order1, order2)
order$dif <- order$M_fun_Z_percentile_1 - order$M_fun_Z_percentile_2
order <- order[order(order$dif),]
dt_plot1$gene_list <- factor(dt_plot1$gene_list, levels = as.character(order$gene_list))

# plot the distribution of relative differnece in funZ for each gene list. 
fig1 <- ggplot(data = dt_plot1, aes(x=gene_list, y=dif, fill = gene_list)) +
  # geom_boxplot() +
  # geom_rect(
  #   aes(xmin = dt$xmin,
  #       xmax = dt$xmax,
  #       ymin = -Inf,
  #       ymax = Inf,
  #       fill = dt$rect_type),
  #   alpha = 0.01) +
  # stat_compare_means(label = "p.signif",
  #                    method = "t.test",
  #                    symnum.args = symnum.args,
#                    size = 6,
#                    hide.ns = T) +
geom_boxplot() +
  # scale_fill_manual(values=c("skyblue")) +
  # scale_fill_brewer(palette="Set2") +
  coord_flip() +
  xlab("") +
  ylab("Differnece in funZ percentile rank\n(mouse orthologue)") +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size = 0.9) +
  annotate("text", x = 1, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 4, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 3, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 2, y = 0.65, label = "****", colour="black", size = 6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14)) 
fig1

##############
### FIGURE ROC
#############


# get predicted data for M_funZ models
sub_score <- subset(dt_long, dt_long$variable == "M_fun_Z_percentile")
terms <- levels(dt_long$gene_list)
Pred.list <- list()

for (i in 1:length(terms)){
  
    sub_term <- subset(sub_score, sub_score$gene_list == terms[i])
    sub_term <- sub_term[!duplicated(sub_term)]
    
    mod <- glm(sig ~ value, data = sub_term, family=binomial)
    
    pred <- predict(mod, sub_term)
    obs <- as.numeric(unlist(sub_term$sig))
    ROC[i] <- auc(roc(obs, pred))
    
    Pred.list[[i]] <- data.frame(H_external_gene_name = sub_term$H_external_gene_name,
                                 gene_list = sub_term$gene_list,
                                 obs = obs,
                                 pred = pred)
}
dt_plot2 <- do.call("rbind", Pred.list)

# add ROC to gene list label
dt_plot2$gene_list <- as.character(dt_plot2$gene_list)
for (cat in unique(dt_plot2$gene_list)){
  ROC <- ROC_table$ROC[ROC_table$gene_list == cat & ROC_table$constraint_score == "M_fun_Z_percentile"]
  dt_plot2$gene_list[dt_plot2$gene_list == cat] <- paste0(dt_plot2$gene_list[dt_plot2$gene_list == cat], " ROC=", ROC)
}
order$gene_list <- as.character(order$gene_list)
for (cat in unique(order$gene_list)){
  ROC <- ROC_table$ROC[ROC_table$gene_list == cat & ROC_table$constraint_score == "M_fun_Z_percentile"]
  order$gene_list[order$gene_list == cat] <- paste0(order$gene_list[order$gene_list == cat], " ROC=", ROC)
}
dt_plot2$gene_list <- factor(dt_plot2$gene_list, levels = as.character(order$gene_list))


fig2 <- ggplot(dt_plot2, aes(d = obs, m = pred, color = gene_list)) + 
  geom_roc(n.cuts = 0) + 
  # scale_colour_brewer(palette="Accent") +
  style_roc() +
  theme(legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.65, 0.25),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
fig2

Fig <- plot_grid(fig1, fig2, ncol = 2, nrow = 1, labels = c("A", "B"), 
                 label_size = 20, hjust = -0.1, vjust = 1.5)
Fig
save_plot("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_gene_list_box_ROC.jpg", Fig, ncol = 1, nrow = 1, base_height = 6, base_width = 13)
fwrite(ROC_table, "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Table_gene_list_ROC.csv")

