### FIGURE 5 -- 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plyr)
library(tidyr)
library(ggpubr)

### IMPORT DATA
data <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")
impc <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
impc.all <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/IMPC_ALL_statistical_results.csv")

### FORMAT

# gene list with no significant phenotype
sig_genes <- unique(impc$marker_symbol)
all_genes <- unique(impc.all$marker_symbol)
no_pheno_genes <- all_genes[-which(all_genes %in% sig_genes)]
no_pheno_genes <- data.frame(external_gene_name = no_pheno_genes,
                             top_level_mp_term_name = rep("no annotated phenotype", length(no_pheno_genes)))


# subset mouse funZ 
data <- data[,c("external_gene_name", "fun_Z_0.0001", "fun_Z_percentile_0.0001")]
# format impc data 
# impc <- impc[impc$zygosity == "homozygote",]
impc <- impc[, c("marker_symbol", "top_level_mp_term_name", "top_level_mp_term_id", "p_value", "effect_size", "percentage_change")]
colnames(impc) <- c("external_gene_name", "top_level_mp_term_name", "top_level_mp_term_id", "p_value", "effect_size", "percentage_change")
# remove genes with no top level phenotype term
impc <- impc[!is.na(impc$top_level_mp_term_name),]
impc <- impc[impc$top_level_mp_term_name != "",]
# split top level phenotype into multiple categories
impc <- separate_rows(impc, top_level_mp_term_name, top_level_mp_term_id, sep = ",")

# Merge
df <- rbind.fill(impc, no_pheno_genes)
df <- merge(df, data)

# remove duplicated top level phenotypes for each gene
temp <- df[,c("external_gene_name", "fun_Z_0.0001", "fun_Z_percentile_0.0001",
              "top_level_mp_term_name", "top_level_mp_term_id")]
temp <- temp[!duplicated(temp),]
# cut top level terms with n < 50
temp <- temp[temp$top_level_mp_term_name %in% names(which(table(temp$top_level_mp_term_name) > 49)), ]
summary(factor(temp$top_level_mp_term_name))

# add reference with all genes
ref <- temp[!duplicated(temp$external_gene_name),]
ref$top_level_mp_term_name <- "all IMPC knockouts"
temp <- rbind.fill(temp, ref)

# calculate funZ percentile for IMPC genes
percentile <- ecdf(temp$fun_Z_0.0001[!duplicated(temp$external_gene_name)])
temp$fun_Z_IMPC_percentile <- percentile(temp$fun_Z_0.0001)
# temp$fun_Z_IMPC_percentile <- ceiling((temp$fun_Z_IMPC_percentile)*100)

# get mean values for all mouse genes
length(data$fun_Z_percentile_0.0001)
y1 <- mean(data$fun_Z_percentile_0.0001, na.rm = T)
# y1 <- median(orths$M_fun_Z_percentile_0.0001, na.rm = T)

# get mean values for all IMPC knockouts
x <- data[data$external_gene_name %in% unique(temp$external_gene_name),]
y2 <- mean(x$fun_Z_percentile_0.0001, na.rm = T)
# y2 <- median(x$M_fun_Z_percentile_0.0001, na.rm = T)

# get mean values for all IMPC knockouts
y3 <- mean(temp$fun_Z_IMPC_percentile[!duplicated(temp$external_gene_name)], na.rm = T)

temp$top_level_mp_term_id[temp$top_level_mp_term_name == "no annotated phenotype"] <- ""
temp$top_level_mp_term_id[temp$top_level_mp_term_name == "all IMPC knockouts"] <- ""
temp <- na.omit(temp)
n.temp <- as.data.frame(table(temp[,"top_level_mp_term_name"]))
colnames(n.temp) <- c("top_level_mp_term_name", "N_genes")
temp <- merge(temp, n.temp)
table(temp$top_level_mp_term_name)
temp$top_level_mp_term_name <- paste0(temp$top_level_mp_term_name, " (n=", temp$N_genes, ")")
# temp$top_level_mp_term_name[temp$top_level_mp_term_name == "no annotated phenotype () (n=1339)"] <- "no annotated phenotype (n=1339)"
# temp$top_level_mp_term_name[temp$top_level_mp_term_name == "all IMPC knockouts () (n=5479)"] <- "all IMPC knockouts (n=5479)"
temp$top_level_mp_term_name <- factor(temp$top_level_mp_term_name)
summary(temp$top_level_mp_term_name)


# set order of top-level terms
# order <- reorder(temp$top_level_mp_term_name, temp$fun_Z_IMPC_percentile, FUN = median, order = T)
order <- aggregate(temp$fun_Z_IMPC_percentile~temp$top_level_mp_term_name, FUN=median)
colnames(order) <- c("top_level_mp_term_name", "fun_Z_IMPC_percentile")
order <- order[order(order$fun_Z_IMPC_percentile),]
order <- order$top_level_mp_term_name[c(2, 1, 3:nrow(order))]
temp$top_level_mp_term_name <- factor(temp$top_level_mp_term_name, levels = order)

# add colour category for baseline
temp$color <- "white"
temp$color[temp$top_level_mp_term_name == "all IMPC knockouts (n=5479)"] <- "blue"

### PLOT IMPC percentile
fig1 <- ggplot(temp, aes(x=top_level_mp_term_name, y=fun_Z_IMPC_percentile, fill = color)) +
  geom_boxplot() +
  scale_fill_manual(values=c("skyblue", "white")) +
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = "all IMPC knockouts (n=5479)",
                     size = 6,
                     hide.ns = T) +
  xlab("Mouse Phenotype Ontology term") +
  ylab("FunZ percentile rank") +
  geom_hline(yintercept=y3, linetype="dashed", color = "red", size = 0.9) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14))
fig1
ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_5AA.jpg", plot = fig1, height = 6, width = 8)
fwrite(temp, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Sandbox/Figure_funZ_IMPC_MP.csv")

### PLOT all mouse gene percentile
fig2 <- ggplot(temp, aes(x=top_level_mp_term_name, y=fun_Z_percentile_0.0001, fill = color)) +
  geom_boxplot() +
  scale_fill_manual(values=c("skyblue", "white")) +
  stat_compare_means(label = "p.signif", method = "wilcox",
                     ref.group = "baseline (n=5479)",
                     size = 6,
                     hide.ns = T) +
  xlab("IMPC top-level phenotype term") +
  ylab("FunZ percentile") +
  geom_hline(yintercept=y1, linetype="dashed", color = "red", size = 0.9) +
  geom_hline(yintercept=y2, linetype="dashed", color = "blue", size = 0.9) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14))
fig2
ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_5B.jpg", plot = fig2, height = 6, width = 8)

### Test difference in funZ between IMPC genes and non-IMPC genes
impc.genes <- unique(temp$external_gene_name)
non.impc.genes <- data$external_gene_name[!data$external_gene_name %in% impc.genes]
wilcox.test(data$fun_Z_percentile_0.0001[data$external_gene_name %in% impc.genes],
            data$fun_Z_percentile_0.0001[data$external_gene_name %in% non.impc.genes])

##############

### Wilcoxon
# test FunZ for phenotype levels against FunZ for all genes, and all genes through IMPC.
terms <- unique(temp$top_level_mp_term_name)
terms <- sort(terms, decreasing = T)
N <- rep(NA, length(terms))
impc.statistic.f <- rep(NA, length(terms))
impc.p.val.f <- rep(NA, length(terms))
means <- rep(NA, length(terms))
medians <- rep(NA, length(terms))

for (i in 1:length(terms)){

  N[i] <- length(temp$fun_Z_IMPC_percentile[temp$top_level_mp_term_name == terms[i]])

  impc.test <- wilcox.test(temp$fun_Z_IMPC_percentile[temp$top_level_mp_term_name == terms[i]], 
                           temp$fun_Z_IMPC_percentile[temp$top_level_mp_term_name == "all IMPC knockouts (n=5479)"])
  impc.statistic.f[i] <- impc.test$statistic
  impc.p.val.f[i] <- impc.test$p.value
  means[i] <- mean(temp$fun_Z_IMPC_percentile[temp$top_level_mp_term_name == terms[i]])
  medians[i] <- median(temp$fun_Z_IMPC_percentile[temp$top_level_mp_term_name == terms[i]])
}
out <- data.frame(MP_term = terms,
                  n_genes = N,
                  mean = means,
                  median = medians,
                  statistic = impc.statistic.f,
                  p_val = impc.p.val.f)

##############


# ### FUNCTIONS
# perc.rank <- function(x) (rank(x)/length(x))*100 
# 
# ## IMPORT DATA
# orths <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
# impc <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
# no.pheno <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/impc_viability_dr7.txt")
# 
# ### FORMAT
# 
# # subset mouse funZ
# orths <- orths[,c("M_external_gene_name", "H_external_gene_name",
#             "M_fun_Z_0.0001", "M_fun_Z_percentile_0.0001")]
# # remove duplicated mouse genes (ie one to many/many to many)
# orths <- orths[!duplicated(orths$M_external_gene_name),]
# 
# # format impc data
# impc <- impc[, c("marker_symbol", "top_level_mp_term_name", "p_value", "effect_size", "percentage_change")]
# colnames(impc) <- c("M_external_gene_name", "top_level_mp_term_name", "p_value", "effect_size", "percentage_change")
# # remove genes with no top level phenotype term
# impc <- impc[!is.na(impc$top_level_mp_term_name),]
# impc <- impc[impc$top_level_mp_term_name != "",]
# # split top level phenotype into multiple categories
# impc <- separate_rows(impc, top_level_mp_term_name, sep = ",")
# 
# # identify genes with no phenotype
# no.pheno <- no.pheno$Approved.symbol[no.pheno$IMPC.Viability.2 == "Viable.Nophenotype"]
# no.pheno <- orths[orths$H_external_gene_name %in% no.pheno,]
# no.pheno <- unique(no.pheno$M_external_gene_name)
# no.pheno <- data.frame(M_external_gene_name = no.pheno,
#                        top_level_mp_term_name = rep("no significant phenotype", length(no.pheno)))
# # hist(no.pheno$M_fun_Z_0.0001)
# # summary(no.pheno$M_fun_Z_0.0001)
# 
# # Merge
# df <- rbind.fill(impc, no.pheno)
# df <- merge(df, orths)
# 
# # remove duplicated top level phenotypes for each gene
# temp <- df[,c("M_external_gene_name", "M_fun_Z_0.0001", "M_fun_Z_percentile_0.0001",
#               "top_level_mp_term_name")]
# temp <- temp[!duplicated(temp),]
# # cut top level terms with n < 50
# temp <- temp[temp$top_level_mp_term_name %in% names(which(table(temp$top_level_mp_term_name) > 49)), ]
# summary(factor(temp$top_level_mp_term_name))
# 
# # add reference with all genes
# ref <- temp[!duplicated(temp$M_external_gene_name),]
# ref$top_level_mp_term_name <- "baseline"
# temp <- rbind.fill(temp, ref)
# 
# ### Wilcoxon
# # test FunZ for phenotype levels against FunZ for all genes, and all genes through IMPC.
# terms <- unique(temp$top_level_mp_term_name)
# N <- rep(NA, length(terms))
# impc.statistic.f <- rep(NA, length(terms))
# impc.p.val.f <- rep(NA, length(terms))
# all.statistic.f <- rep(NA, length(terms))
# all.p.val.f <- rep(NA, length(terms))
# 
# for (i in 1:length(terms)){
# 
#   N[i] <- length(temp$M_fun_Z_percentile_0.0001[temp$top_level_mp_term_name == terms[i]])
# 
#   impc.test <- wilcox.test(temp$M_fun_Z_percentile_0.0001[temp$top_level_mp_term_name == terms[i]], orths$M_fun_Z_percentile_0.0001[orths$M_external_gene_name %in% unique(temp$M_external_gene_name)])
#   impc.statistic.f[i] <- impc.test$statistic
#   impc.p.val.f[i] <- impc.test$p.value
# 
#   all.test <- wilcox.test(temp$M_fun_Z_percentile_0.0001[temp$top_level_mp_term_name == terms[i]], orths$M_fun_Z_percentile_0.0001)
#   all.statistic.f[i] <- all.test$statistic
#   all.p.val.f[i] <- all.test$p.value
# }
# out <- data.frame(Term = terms,
#                   n_genes = N,
#                   statistic_all_FunZ = all.statistic.f,
#                   p_all_FunZ = all.p.val.f,
#                   statistic_impc_FunZ = impc.statistic.f,
#                   p_impc_FunZ = impc.p.val.f)
# 
# # get mean values for all mouse genes
# y1 <- mean(orths$M_fun_Z_percentile_0.0001, na.rm = T)
# # y1 <- median(orths$M_fun_Z_percentile_0.0001, na.rm = T)
# 
# # get mean values for all IMPC knockouts
# x <- orths[orths$M_external_gene_name %in% unique(temp$M_external_gene_name),]
# y2 <- mean(x$M_fun_Z_percentile_0.0001, na.rm = T)
# # y2 <- median(x$M_fun_Z_percentile_0.0001, na.rm = T)
# 
# temp <- na.omit(temp)
# n.temp <- as.data.frame(table(temp[,"top_level_mp_term_name"]))
# colnames(n.temp) <- c("top_level_mp_term_name", "N_genes")
# temp <- merge(temp, n.temp)
# temp$top_level_mp_term_name <- paste0(temp$top_level_mp_term_name, " (n=", temp$N_genes, ")")
# summary(factor(temp$top_level_mp_term_name))
# 
# fig <- ggplot(temp, aes(x=reorder(top_level_mp_term_name, M_fun_Z_percentile_0.0001, FUN = median, order = T), y=M_fun_Z_percentile_0.0001)) +
#   geom_boxplot() +
#   stat_compare_means(label = "p.signif", method = "wilcox",
#                      ref.group = "baseline (n=3750)",
#                      size = 6,
#                      hide.ns = T) +
#   xlab("IMPC top-level phenotype term") +
#   ylab("FunZ percentile") +
#   geom_hline(yintercept=y1, linetype="dashed", color = "red") +
#   geom_hline(yintercept=y2, linetype="dashed", color = "blue") +
#   coord_flip() +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         text = element_text(size = 14))
# fig
# ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_5.jpg", plot = fig, height = 6, width = 8)
# 
# #####
# 
# # ### IMPORT
# # impc <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/ALL_genotype_phenotype.csv")
# # newC <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")
# # 
# # ### FORMAT
# # 
# # # colnames
# # 
# # 
# # # percentiles
# # RVIS_esp$RVIS_ESP_percentile <- RVIS_esp$RVIS_ESP_percentile/100
# # percentile <- ecdf(Lek$mis_z_Lek)
# # Lek$mis_z_percentile_Lek <- percentile(Lek$mis_z_Lek)
# # percentile <- ecdf(Lek$pLI_Lek)
# # Lek$pLI_percentile_Lek <- percentile(Lek$pLI_Lek)
# # 
# # df <- merge(Lek, RVIS_esp, all = T)
# # df <- merge(df, newC, all = T)
# # df <- df[!is.na(df$H_external_gene_name),]
# # df <- df[df$H_external_gene_name != "",]
# # 
# # impc <- df[df$H_external_gene_name %in% no.pheno$Approved.symbol,]
# # impc$cat <- "Phenotype"
# # impc$cat[impc$H_external_gene_name %in% no.pheno$Approved.symbol[no.pheno$IMPC.Viability.2 == "Lethal"]] <- "Lethal"
# # impc$cat[impc$H_external_gene_name %in% no.pheno$Approved.symbol[no.pheno$IMPC.Viability.2 == "Viable.Nophenotype"]] <- "No_phenotype"
# # table(impc$cat)
# # 
# # impc$Lethal <- 0
# # impc$Lethal[impc$cat == "Lethal"] <- 1
# # glm <- glm(Lethal ~ RVIS_ESP_percentile, data = impc, family=binomial)
# # summary(glm)
# # pred <- predict(glm, impc)
# # obs <- impc$Lethal
# # auc(roc(obs, pred))
# 

# out <- temp[,c("top_level_mp_term_name",  "external_gene_name", "fun_Z_IMPC_percentile")]
# library(tidyr)
# library(reshape2)
# data_wide <- dcast(out, external_gene_name~top_level_mp_term_name, value.var="fun_Z_IMPC_percentile")
# fwrite(data_wide, "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Sandbox/IMPC_plot_wide.csv")
