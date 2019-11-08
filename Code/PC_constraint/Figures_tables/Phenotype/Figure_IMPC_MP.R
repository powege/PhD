rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr) 
library(ggpubr)

### IMPORT
# procedures and life stage
PPI <- fread("~/Dropbox/PhD/Data/IMPC/procedure_pipeline_info.csv")
# procedures for each knockout
LPI <- fread("~/Dropbox/PhD/Data/IMPC/Line_procedure_info.csv")
# significant MP terms
impc <- fread("~/Dropbox/PhD/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
# all statistical results
impc.all <- fread("~/Dropbox/PhD/Data/IMPC/IMPC_ALL_statistical_results.csv")
# funZ
funz <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv")

# create long dt knockout; procedure key; life stage
# sub <- subset (LPI, LPI$gene_symbol == "Xbp1")
wolf <- function(sub){
  procedure_keys <- unlist(paste(sub$procedure_keys, collapse = ","))
  procedure_keys <- unique(unlist(strsplit(procedure_keys, ",")))
  out <- data.table(procedure_key = procedure_keys)
  return(out)
}
gene_procedures <- ddply(LPI, "gene_symbol", wolf)
lifestage <- PPI[,c("procedure_key", "life_stage")]
gene_procedures <- lifestage[gene_procedures, on = "procedure_key"]
colnames(gene_procedures) <- c("procedure_stable_id", "life_stage", "gene_symbol" )
gene_procedures <- gene_procedures[!duplicated(gene_procedures),]

# create dt of all potential top-level MP for each gene
dt.all <- impc.all[,c("marker_symbol", "procedure_stable_id", "top_level_mp_term_name")]
colnames(dt.all) <- c("gene_symbol", "procedure_stable_id", "top_level_mp_term_name")
# remove genes with no top level phenotype term
dt.all <- dt.all[!is.na(dt.all$top_level_mp_term_name),]
dt.all <- dt.all[dt.all$top_level_mp_term_name != "",]
# split top level phenotype into multiple categories
dt.all <- separate_rows(dt.all, top_level_mp_term_name, sep = ",")
dt.all <- dt.all[!duplicated(dt.all),]

# create dt  of significant top-level MP for each gene
dt.sig <- impc[,c("marker_symbol", "procedure_stable_id", "top_level_mp_term_name")]
colnames(dt.sig) <- c("gene_symbol", "procedure_stable_id", "top_level_mp_term_name")
# remove genes with no top level phenotype term
dt.sig <- dt.sig[!is.na(dt.sig$top_level_mp_term_name),]
dt.sig <- dt.sig[dt.sig$top_level_mp_term_name != "",]
# split top level phenotype into multiple categories
dt.sig <- separate_rows(dt.sig, top_level_mp_term_name, sep = ",")
dt.sig <- dt.sig[!duplicated(dt.sig),]

# identify all potential MP terms per knockout
dt.all <- dt.all[,c("top_level_mp_term_name", "gene_symbol")]
dt.all <- dt.all[!duplicated(dt.all),]
# cf.all <- subset(dt.all, dt.all$top_level_mp_term_name == "embryo phenotype")

# identify significant MP temrs per knockout
dt.sig <- dt.sig[,c("top_level_mp_term_name", "gene_symbol")]
dt.sig <- dt.sig[!duplicated(dt.sig),]
# cf.sig <- subset(dt.sig, dt.sig$top_level_mp_term_name == "embryo phenotype")

# ientify non-significant MP terms per knockout 
dt.nosig <- anti_join(dt.all, dt.sig)
# nrow(dt.all) - nrow(dt.sig)

# add significance and merge
dt.sig$sig <- 1
dt.nosig$sig <- 0
dt <- rbind(dt.sig, dt.nosig)
# t1 <- paste0(dt.sig$top_level_mp_term_name, dt.sig$gene_symbol)
# t2 <- paste0(dt.nosig$top_level_mp_term_name, dt.nosig$gene_symbol)
# summary(t2 %in% t1)

# add gene list with no significant phenotype
sig_genes <- unique(impc$marker_symbol)
all_genes <- unique(impc.all$marker_symbol)
no_pheno_genes <- all_genes[-which(all_genes %in% sig_genes)]
no_pheno_genes <- data.frame(gene_symbol = c(no_pheno_genes, sig_genes),
                             top_level_mp_term_name = c(rep("no annotated phenotype", length(no_pheno_genes)), rep("no annotated phenotype", length(sig_genes))),
                             sig = c(rep(1, length(no_pheno_genes)), rep(0, length(sig_genes))))
dt <- rbind(dt, no_pheno_genes)

# format funZ 
funz <- funz[,c("external_gene_name", "fun_Z")]
# funz <- funz[,c("external_gene_name", "fun_Z_MSE")]
colnames(funz) <- c("gene_symbol", "funZ")

# calculate funZ percentile for all mouse genes
percentile <- ecdf(funz$funZ[!duplicated(funz$funZ)])
funz$funZ_IMPC_percentile <- percentile(funz$funZ)

# merge with funZ 
dt <- merge(dt, funz)
length(unique(dt$gene_symbol))
# summary(unique(dt$gene_symbol) %in% unique(funz$gene_symbol))

# # calculate funZ percentile for IMPC genes
# percentile <- ecdf(dt$funZ[!duplicated(dt$funZ)])
# dt$funZ_IMPC_percentile <- percentile(dt$funZ)


# cut top level terms with n < 50
temp <- dt[!duplicated(dt),]
tmp <- subset(temp, temp$sig == 1)
temp <- temp[temp$top_level_mp_term_name %in% names(which(table(tmp$top_level_mp_term_name) > 49)), ]
N_sig <- table(temp[,c("top_level_mp_term_name", "sig")])

# plot formatting 
# set order of top-level terms
tmp1 <- subset(temp, temp$sig == 1)
order1 <- aggregate(tmp1$funZ_IMPC_percentile~tmp1$top_level_mp_term_name, FUN=median)
tmp2 <- subset(temp, temp$sig == 0)
order2 <- aggregate(tmp2$funZ_IMPC_percentile~tmp2$top_level_mp_term_name, FUN=median)
colnames(order1) <- c("top_level_mp_term_name", "funZ_IMPC_percentile_1")
colnames(order2) <- c("top_level_mp_term_name", "funZ_IMPC_percentile_2")
order <- merge(order1, order2)
order$dif <- order$funZ_IMPC_percentile_1 - order$funZ_IMPC_percentile_2
order <- order[order(order$dif),]
temp$top_level_mp_term_name <- factor(temp$top_level_mp_term_name, levels = as.character(order$top_level_mp_term_name))
temp$sig <- as.factor(temp$sig)
y.int <- median(temp$funZ_IMPC_percentile[!duplicated(temp$gene_symbol)], na.rm = T)

# calculate median funZ percentile for each KO null group
null_average <- subset(temp, temp$sig == 0)
null_average <- aggregate(funZ_IMPC_percentile ~ top_level_mp_term_name, data=null_average, FUN=median)
colnames(null_average) <- c("top_level_mp_term_name", "null_average")

# calculate the relative difference in funZ for each KO group relative to null average
plot_df <- subset(temp, temp$sig == 1)
plot_df <- merge(plot_df, null_average, by = "top_level_mp_term_name")
plot_df$dif <- plot_df$funZ_IMPC_percentile - plot_df$null_average


# Mann Whitney U test differneces in funZ 
terms <- levels(temp$top_level_mp_term_name)
N_sig <- rep(NA, length(terms))
N_nosig <- rep(NA, length(terms))
impc.statistic.f <- rep(NA, length(terms))
impc.p.val.f <- rep(NA, length(terms))
median_sig <- rep(NA, length(terms))
median_nosig <- rep(NA, length(terms))

for (i in 1:length(terms)){
  
  N_sig[i] <- length(temp$gene_symbol[temp$top_level_mp_term_name == terms[i] & temp$sig == 1])
  N_nosig[i] <- length(temp$gene_symbol[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
  
  impc.test <- wilcox.test(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 1], 
                           temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
  impc.statistic.f[i] <- impc.test$statistic
  impc.p.val.f[i] <- impc.test$p.value
  median_sig[i] <- median(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 1])
  median_nosig[i] <- median(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
}
MWU_table <- data.frame(
  top_level_MP_term = terms,
  n_KO_with_MP = N_sig,
  n_KO_without_MP = N_nosig,
  median_funZ_with_MP = round(median_sig, digits = 3),
  median_funZ_without_MP = round(median_nosig, digits = 3),
  statistic = formatC(impc.statistic.f, format = "e", digits = 2),
  p_val = formatC(p.adjust(impc.p.val.f, method = "bonferroni", n = length(impc.p.val.f)), format = "e", digits = 2)
)
MWU_table$funZ_diff <- MWU_table$median_funZ_with_MP - MWU_table$median_funZ_without_MP



# plot the distribution of relative differnece in funZ for each KO group. 
fig <- ggplot(data = plot_df, aes(x=top_level_mp_term_name, y=dif)) +
  # geom_boxplot() +
  # geom_rect(
  #   aes(xmin = temp$xmin,
  #       xmax = temp$xmax,
  #       ymin = -Inf,
  #       ymax = Inf,
  #       fill = temp$rect_type),
  #   alpha = 0.01) +
  # stat_compare_means(label = "p.signif",
  #                    method = "t.test",
  #                    symnum.args = symnum.args,
  #                    size = 6,
  #                    hide.ns = T) +
  geom_boxplot(fill = "skyblue") +
  # scale_fill_manual(values=c("skyblue")) +
  coord_flip() +
  xlab("Mouse Phenotype Ontology term") +
  ylab("Differnece in funZ percentile rank") +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size = 0.9) +
  annotate("text", x = 22, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 21, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 20, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 18, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 16, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 14, y = 0.65, label = "**", colour="black", size = 6) +
  annotate("text", x = 13, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 12, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 11, y = 0.65, label = "**", colour="black", size = 6) +
  annotate("text", x = 1, y = 0.65, label = "****", colour="black", size = 6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14)) 
fig


### OUTPUT
ggsave("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_IMPC_MP.jpg", plot = fig, height = 6.5, width = 7.5)
fwrite(MWU_table, "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Table_IMPC_MP.csv")



#####

# ### RUN WITH DIFFERNT CONSTRAINT METRIC
# 
# ### IMPORT
# # procedures and life stage
# PPI <- fread("~/Dropbox/PhD/Data/IMPC/procedure_pipeline_info.csv")
# # procedures for each knockout
# LPI <- fread("~/Dropbox/PhD/Data/IMPC/Line_procedure_info.csv")
# # significant MP terms
# impc <- fread("~/Dropbox/PhD/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
# # all statistical results
# impc.all <- fread("~/Dropbox/PhD/Data/IMPC/IMPC_ALL_statistical_results.csv")
# # funZ
# funz <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv")
# 
# # create long dt knockout; procedure key; life stage
# # sub <- subset (LPI, LPI$gene_symbol == "Xbp1")
# wolf <- function(sub){
#   procedure_keys <- unlist(paste(sub$procedure_keys, collapse = ","))
#   procedure_keys <- unique(unlist(strsplit(procedure_keys, ",")))
#   out <- data.table(procedure_key = procedure_keys)
#   return(out)
# }
# gene_procedures <- ddply(LPI, "gene_symbol", wolf)
# lifestage <- PPI[,c("procedure_key", "life_stage")]
# gene_procedures <- lifestage[gene_procedures, on = "procedure_key"]
# colnames(gene_procedures) <- c("procedure_stable_id", "life_stage", "gene_symbol" )
# gene_procedures <- gene_procedures[!duplicated(gene_procedures),]
# 
# # create dt of all potential top-level MP for each gene
# dt.all <- impc.all[,c("marker_symbol", "procedure_stable_id", "top_level_mp_term_name")]
# colnames(dt.all) <- c("gene_symbol", "procedure_stable_id", "top_level_mp_term_name")
# # remove genes with no top level phenotype term
# dt.all <- dt.all[!is.na(dt.all$top_level_mp_term_name),]
# dt.all <- dt.all[dt.all$top_level_mp_term_name != "",]
# # split top level phenotype into multiple categories
# dt.all <- separate_rows(dt.all, top_level_mp_term_name, sep = ",")
# dt.all <- dt.all[!duplicated(dt.all),]
# 
# # create dt  of significant top-level MP for each gene
# dt.sig <- impc[,c("marker_symbol", "procedure_stable_id", "top_level_mp_term_name")]
# colnames(dt.sig) <- c("gene_symbol", "procedure_stable_id", "top_level_mp_term_name")
# # remove genes with no top level phenotype term
# dt.sig <- dt.sig[!is.na(dt.sig$top_level_mp_term_name),]
# dt.sig <- dt.sig[dt.sig$top_level_mp_term_name != "",]
# # split top level phenotype into multiple categories
# dt.sig <- separate_rows(dt.sig, top_level_mp_term_name, sep = ",")
# dt.sig <- dt.sig[!duplicated(dt.sig),]
# 
# # identify all potential MP terms per knockout
# dt.all <- dt.all[,c("top_level_mp_term_name", "gene_symbol")]
# dt.all <- dt.all[!duplicated(dt.all),]
# # cf.all <- subset(dt.all, dt.all$top_level_mp_term_name == "embryo phenotype")
# 
# # identify significant MP temrs per knockout
# dt.sig <- dt.sig[,c("top_level_mp_term_name", "gene_symbol")]
# dt.sig <- dt.sig[!duplicated(dt.sig),]
# # cf.sig <- subset(dt.sig, dt.sig$top_level_mp_term_name == "embryo phenotype")
# 
# # ientify non-significant MP terms per knockout
# dt.nosig <- anti_join(dt.all, dt.sig)
# # nrow(dt.all) - nrow(dt.sig)
# 
# # add significance and merge
# dt.sig$sig <- 1
# dt.nosig$sig <- 0
# dt <- rbind(dt.sig, dt.nosig)
# # t1 <- paste0(dt.sig$top_level_mp_term_name, dt.sig$gene_symbol)
# # t2 <- paste0(dt.nosig$top_level_mp_term_name, dt.nosig$gene_symbol)
# # summary(t2 %in% t1)
# 
# # add gene list with no significant phenotype
# sig_genes <- unique(impc$marker_symbol)
# all_genes <- unique(impc.all$marker_symbol)
# no_pheno_genes <- all_genes[-which(all_genes %in% sig_genes)]
# no_pheno_genes <- data.frame(gene_symbol = c(no_pheno_genes, sig_genes),
#                              top_level_mp_term_name = c(rep("no annotated phenotype", length(no_pheno_genes)), rep("no annotated phenotype", length(sig_genes))),
#                              sig = c(rep(1, length(no_pheno_genes)), rep(0, length(sig_genes))))
# dt <- rbind(dt, no_pheno_genes)
# 
# # format funZ
# funz <- funz[,c("external_gene_name", "fun_Z_adj")]
# colnames(funz) <- c("gene_symbol", "funZ")
# 
# # calculate funZ percentile for all mouse genes
# percentile <- ecdf(funz$funZ[!duplicated(funz$funZ)])
# funz$funZ_IMPC_percentile <- percentile(funz$funZ)
# 
# # merge with funZ
# dt <- merge(dt, funz)
# length(unique(dt$gene_symbol))
# # summary(unique(dt$gene_symbol) %in% unique(funz$gene_symbol))
# 
# # # calculate funZ percentile for IMPC genes
# # percentile <- ecdf(dt$funZ[!duplicated(dt$funZ)])
# # dt$funZ_IMPC_percentile <- percentile(dt$funZ)
# 
# 
# # cut top level terms with n < 50
# temp <- dt[!duplicated(dt),]
# tmp <- subset(temp, temp$sig == 1)
# temp <- temp[temp$top_level_mp_term_name %in% names(which(table(tmp$top_level_mp_term_name) > 49)), ]
# N_sig <- table(temp[,c("top_level_mp_term_name", "sig")])
# 
# # plot formatting
# # set order of top-level terms
# tmp1 <- subset(temp, temp$sig == 1)
# order1 <- aggregate(tmp1$funZ_IMPC_percentile~tmp1$top_level_mp_term_name, FUN=median)
# tmp2 <- subset(temp, temp$sig == 0)
# order2 <- aggregate(tmp2$funZ_IMPC_percentile~tmp2$top_level_mp_term_name, FUN=median)
# colnames(order1) <- c("top_level_mp_term_name", "funZ_IMPC_percentile_1")
# colnames(order2) <- c("top_level_mp_term_name", "funZ_IMPC_percentile_2")
# order <- merge(order1, order2)
# order$dif <- order$funZ_IMPC_percentile_1 - order$funZ_IMPC_percentile_2
# order <- order[order(order$dif),]
# temp$top_level_mp_term_name <- factor(temp$top_level_mp_term_name, levels = as.character(order$top_level_mp_term_name))
# temp$sig <- as.factor(temp$sig)
# y.int <- median(temp$funZ_IMPC_percentile[!duplicated(temp$gene_symbol)], na.rm = T)
# 
# # calculate median funZ percentile for each KO null group
# null_average <- subset(temp, temp$sig == 0)
# null_average <- aggregate(funZ_IMPC_percentile ~ top_level_mp_term_name, data=null_average, FUN=median)
# colnames(null_average) <- c("top_level_mp_term_name", "null_average")
# 
# # calculate the relative difference in funZ for each KO group relative to null average
# plot_df <- subset(temp, temp$sig == 1)
# plot_df <- merge(plot_df, null_average, by = "top_level_mp_term_name")
# plot_df$dif <- plot_df$funZ_IMPC_percentile - plot_df$null_average
# 
# 
# # Mann Whitney U test differneces in funZ
# terms <- levels(temp$top_level_mp_term_name)
# N_sig <- rep(NA, length(terms))
# N_nosig <- rep(NA, length(terms))
# impc.statistic.f <- rep(NA, length(terms))
# impc.p.val.f <- rep(NA, length(terms))
# median_sig <- rep(NA, length(terms))
# median_nosig <- rep(NA, length(terms))
# 
# for (i in 1:length(terms)){
# 
#   N_sig[i] <- length(temp$gene_symbol[temp$top_level_mp_term_name == terms[i] & temp$sig == 1])
#   N_nosig[i] <- length(temp$gene_symbol[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
# 
#   impc.test <- wilcox.test(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 1],
#                            temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
#   impc.statistic.f[i] <- impc.test$statistic
#   impc.p.val.f[i] <- impc.test$p.value
#   median_sig[i] <- median(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 1])
#   median_nosig[i] <- median(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
# }
# MWU_table <- data.frame(
#   top_level_MP_term = terms,
#   n_KO_with_MP = N_sig,
#   n_KO_without_MP = N_nosig,
#   median_funZ_with_MP = round(median_sig, digits = 3),
#   median_funZ_without_MP = round(median_nosig, digits = 3),
#   statistic = formatC(impc.statistic.f, format = "e", digits = 2),
#   p_val = formatC(p.adjust(impc.p.val.f, method = "bonferroni", n = length(impc.p.val.f)), format = "e", digits = 2)
# )
# MWU_table$funZ_diff <- MWU_table$median_funZ_with_MP - MWU_table$median_funZ_without_MP
# 
# 
# 
# # plot the distribution of relative differnece in funZ for each KO group.
# fig <- ggplot(data = plot_df, aes(x=top_level_mp_term_name, y=dif)) +
#   # geom_boxplot() +
#   # geom_rect(
#   #   aes(xmin = temp$xmin,
#   #       xmax = temp$xmax,
#   #       ymin = -Inf,
#   #       ymax = Inf,
#   #       fill = temp$rect_type),
#   #   alpha = 0.01) +
#   # stat_compare_means(label = "p.signif",
#   #                    method = "t.test",
#   #                    symnum.args = symnum.args,
# #                    size = 6,
# #                    hide.ns = T) +
# geom_boxplot(fill = "skyblue") +
#   # scale_fill_manual(values=c("skyblue")) +
#   coord_flip() +
#   xlab("Mouse Phenotype Ontology term") +
#   ylab("Differnece in funZ percentile rank") +
#   geom_hline(yintercept=0, linetype="dashed", color = "red", size = 0.9) +
#   # annotate("text", x = 22, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 21, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 20, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 18, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 16, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 14, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 13, y = 0.65, label = "***", colour="black", size = 6) +
#   # annotate("text", x = 12, y = 0.65, label = "**", colour="black", size = 6) +
#   # annotate("text", x = 11, y = 0.65, label = "****", colour="black", size = 6) +
#   # annotate("text", x = 1, y = 0.65, label = "****", colour="black", size = 6) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "none",
#         text = element_text(size = 14))
# fig


### OUTPUT
# ggsave("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_IMPC_MP.jpg", plot = fig, height = 6.5, width = 7.5)
# fwrite(MWU_table, "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Table_IMPC_MP.csv")


