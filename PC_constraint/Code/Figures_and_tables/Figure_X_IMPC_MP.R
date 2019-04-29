rm(list = ls())
graphics.off()

library(data.table)
library(plyr)
library(dplyr) 
library(ggpubr)

### IMPORT
# procedures and life stage
PPI <- fread("~/Downloads/procedure_pipeline_info.csv")
# procedures for each knockout
LPI <- fread("~/Downloads/Line_procedure_info.csv")
# significant MP terms
impc <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
# all statistical results
impc.all <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/IMPC_ALL_statistical_results.csv")
# funZ
funz <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")

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

# merge with funZ 
funz <- funz[,c("external_gene_name", "fun_Z_0.0001")]
colnames(funz) <- c("gene_symbol", "funZ")
dt <- merge(dt, funz)
# length(unique(dt$gene_symbol))
# summary(unique(dt$gene_symbol) %in% unique(funz$gene_symbol))

# calculate funZ percentile for IMPC genes
percentile <- ecdf(dt$funZ[!duplicated(dt$funZ)])
dt$funZ_IMPC_percentile <- percentile(dt$funZ)
# hist(dt$funZ_IMPC_percentile)

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
y.int <- mean(temp$funZ_IMPC_percentile[!duplicated(temp$gene_symbol)], na.rm = T)


# Bonferoni correction
N = 22 # Number of comparisons
symnum.args <- list(
  cutpoints = c(0, 0.0001/N, 0.001/N, 0.01/N, 0.05/N, 1),
  symbols = c("****", "***", "**", "*", "ns")
  # symbols = c("<0.0001", "<0.001", "<0.01", "<0.05", "ns")
)

rects <- data.frame(xmin = head(seq <- seq(0.5, 22 + .5, 1), -1), 
                    xmax = tail(seq, -1), 
                    top_level_mp_term_name = levels(temp$top_level_mp_term_name),
                    rect_type = c("a", "c"))
temp <- merge(temp, rects, by = "top_level_mp_term_name")

# plot 
fig <- ggplot(data = temp, aes(x=top_level_mp_term_name, y=funZ_IMPC_percentile, fill = sig)) +
  geom_boxplot() +
  geom_rect(
            aes(xmin = temp$xmin,
                xmax = temp$xmax,
                ymin = -Inf,
                ymax = Inf,
                fill = temp$rect_type),
            alpha = 0.01) +
  stat_compare_means(label = "p.signif", 
                     method = "wilcox",
                     symnum.args = symnum.args, 
                     size = 6,
                     hide.ns = T) +
  geom_boxplot() +
  scale_fill_manual(values=c("skyblue", "white", "grey", "white")) +
  coord_flip() +
  xlab("Mouse Phenotype Ontology term") +
  ylab("FunZ percentile rank") +
  geom_hline(yintercept=y.int, linetype="dashed", color = "red", size = 0.9) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14)) 
fig
ggsave("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/Figure_X_IMPC_MP.jpg", plot = fig, height = 8, width = 8)
length(unique(temp$gene_symbol)) # N genes under analysis
length(unique(temp$gene_symbol[temp$top_level_mp_term_name == "no annotated phenotype" & temp$sig == 1]))

##################

# calculate median funZ percentile for each KO null group
null_average <- subset(temp, temp$sig == 0)
null_average <- aggregate(funZ_IMPC_percentile ~ top_level_mp_term_name, data=null_average, FUN=median)
colnames(null_average) <- c("top_level_mp_term_name", "null_average")

# calculate the relative difference in funZ for each KO group relative to null median
plot_df <- subset(temp, temp$sig == 1)
plot_df <- merge(plot_df, null_average)
plot_df$dif <- plot_df$funZ_IMPC_percentile - plot_df$null_average

# plot the distribution of relative differnece in funZ for each KO group. 
# ADD SIGNIFICANCE
fig2 <- ggplot(data = plot_df, aes(x=top_level_mp_term_name, y=dif)) +
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
  annotate("text", x = 19, y = 0.65, label = "**", colour="black", size = 6) +
  annotate("text", x = 18, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 17, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 16, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 14, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 13, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 11, y = 0.65, label = "**", colour="black", size = 6) +
  annotate("text", x = 10, y = 0.65, label = "****", colour="black", size = 6) +
  annotate("text", x = 1, y = 0.65, label = "****", colour="black", size = 6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 14)) 
fig2
ggsave("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/Figure_X_IMPC_MP.jpg", plot = fig2, height = 6.5, width = 7.5)



# check number of procedures per lifestage for each top-level MP
# tmp <- dt.all[,c("procedure_stable_id","top_level_mp_term_name", "life_stage")]
# tmp <- tmp[!duplicated(tmp),]
# tmp <- tmp[,c("top_level_mp_term_name", "life_stage")]
# tmp <- as.data.table(table(tmp))


