### Relationship between constraint and number of pheotypes (total, and top-level, with prenatal phenotypes removed)
### Test relationship in contraint between genes that cause prenatal lethality and preweaning lethality
### Test for differneces in constraint between top-level phenotype groups. 
### Data from IMPC release 9.1 (http://www.mousephenotype.org/data/release)

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(tidyr)
library("cowplot")

### ROMOVE GENES WITH PRENATAL LETHAL PHENOTYPE
### test for differnece in funZ between genes that cause pre weaning lethality and those that dont
### Plot mean/median counts for all and top only phenotypes

### FUNCTION for plotting lm equation on ggplot
linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 2),
            yy = format(abs(coef(k)[2]), digits = 2),
            r2 = format(summary(k)$r.squared, digits = 3));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)   
  }
  as.character(as.expression(eq));               
}

### IMPORT DATA
data <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")
impc <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
impc.all <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/IMPC_ALL_statistical_results.csv")

### FORMAT

# gene list of prenatal and embryoniic lethal
lethal_genes <- impc[grep("lethal", impc$mp_term_name),]
# table(lethal_genes$mp_term_name)
embryonic_lethal <- unique(impc$marker_symbol[grep("embryonic lethal", impc$mp_term_name)])
prenatal_lethal <- unique(impc$marker_symbol[grep("prenatal lethal", impc$mp_term_name)])
prenatal_genes <- unique(c(embryonic_lethal, prenatal_lethal))

# gene list with no significant phenotype
sig_genes <- unique(impc$marker_symbol)
all_genes <- unique(impc.all$marker_symbol)
no_pheno_genes <- all_genes[-which(all_genes %in% sig_genes)]
no_pheno_genes <- data.frame(external_gene_name = no_pheno_genes,
                             pheno_count_all = rep(0, length(no_pheno_genes)),
                             pheno_count_top = rep(0, length(no_pheno_genes)))


# subset mouse funZ 
data <- data[,c("external_gene_name", "fun_Z_0.0001", "fun_Z_percentile_0.0001")]
# format impc data 
impc <- impc[, c("marker_symbol", "top_level_mp_term_name", "p_value", "effect_size", "percentage_change")]
colnames(impc) <- c("external_gene_name", "top_level_mp_term_name", "p_value", "effect_size", "percentage_change")
# remove genes with no top level phenotype term
impc <- impc[!is.na(impc$top_level_mp_term_name),]
impc <- impc[impc$top_level_mp_term_name != "",]
# split top level phenotype into multiple categories
impc <- separate_rows(impc, top_level_mp_term_name, sep = ",")
# rempve % from percentage_cjhange
impc$percentage_change = substr(impc$percentage_change, 1, nchar(impc$percentage_change)-1)
impc$percentage_change <- as.numeric(impc$percentage_change)

# count top-level terms per gene
all.counts <- as.data.frame(table(impc$external_gene_name))
colnames(all.counts) <- c("external_gene_name", "pheno_count_all")

# remove duplicated top level phenotypes for each gene
impc_rem <- impc[!duplicated(impc[,c("external_gene_name","top_level_mp_term_name")]),]
# count top-level terms per gene
top.counts <- as.data.frame(table(impc_rem$external_gene_name))
colnames(top.counts) <- c("external_gene_name", "pheno_count_top")

# Merge, add no pheno genes, remove prenatal lethal genes
df <- merge(top.counts, all.counts)
df <- df[-which(df$external_gene_name %in% prenatal_genes),]
df <- rbind(df, no_pheno_genes)
df <- merge(df, data)

# calculate funZ percentile for IMPC genes
percentile <- ecdf(df$fun_Z_0.0001)
df$fun_Z_IMPC_percentile <- percentile(df$fun_Z_0.0001)
df$fun_Z_IMPC_percentile <- ceiling((df$fun_Z_IMPC_percentile)*100)
# hist(df$fun_Z_IMPC_percentile)
# hist(df$fun_Z_percentile_0.0001)

# count mean and median phenotypes per percentile
df_mn_all <- aggregate(df$pheno_count_all, by=list(Category=df$fun_Z_IMPC_percentile), FUN=mean)
df_md_all <- aggregate(df$pheno_count_all, by=list(Category=df$fun_Z_IMPC_percentile), FUN=median)
df_mn_top <- aggregate(df$pheno_count_top, by=list(Category=df$fun_Z_IMPC_percentile), FUN=mean)
df_md_top <- aggregate(df$pheno_count_top, by=list(Category=df$fun_Z_IMPC_percentile), FUN=median)

# colnames(df_mn_all) <- c("FunZ_percentile", "n_pheno")
# colnames(df_mn_top) <- c("FunZ_percentile", "n_top_pheno")
# df_plot <- merge(df_mn_all, df_mn_top)
# df_plot <- melt(df_plot, id.vars = "FunZ_percentile")
  
  mod_all <- lm(df_md_all$x~df_md_all$Category)
  equation_all <- linear(mod_all)

  plot_all <- ggplot(df_md_all, aes(x = Category, y = x)) +
    geom_point() +
    geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
    annotate("text", x = 37, y = 7.5, label = equation_all, colour="black", size = 5, parse=TRUE) +
    xlab("FunZ percentile") +
    ylab('Median number of phenotypes') +
    # xlim(-10, 10) +
    # ylim(-10, 10) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 14),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          plot.margin=unit(c(1.5,1,1,1),"cm"))
  plot_all
  
  mod_top <- lm(df_md_top$x~df_md_top$Category)
  equation_top <- linear(mod_top)
  
  plot_top <- ggplot(df_md_top, aes(x = Category, y = x)) +
    geom_point() +
    geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
    annotate("text", x = 40, y = 3.1, label = equation_top, colour="black", size = 5, parse=TRUE) +
    xlab("FunZ percentile") +
    ylab('Median number of top-level phenotypes') +
    # xlim(-10, 10) +
    # ylim(-10, 10) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 14),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          plot.margin=unit(c(1.5,1,1,1),"cm"))
  plot_top

  Fig <- plot_grid(plot_all, plot_top, ncol = 2, nrow = 1, labels = c("A", "B"), 
                   label_size = 20, hjust = -0.1, vjust = 1.5)
  Fig
  save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_pleiotropy.jpg", Fig, ncol = 1, nrow = 1, base_height = 5, base_width = 10)
  
  
  
#######

# df2 <- merge(impc, data, all = T)
# df2 <- na.omit(df2)
# 
# percentile <- ecdf(df2$fun_Z_0.0001)
# df2$fun_Z_IMPC_percentile <- percentile(df2$fun_Z_0.0001)
# df2$fun_Z_IMPC_percentile <- ceiling((df2$fun_Z_IMPC_percentile)*100)
# 
# 
# df2 <- subset(df2, df2$percentage_change > -1000)
# 
# df2$percentage_change_abs <- abs(df2$percentage_change)
# 
# plot(df2$percentage_change_abs ~ df2$fun_Z_IMPC_percentile)
# mod <- lm(df2$percentage_change_abs ~ df2$fun_Z_IMPC_percentile + df2$top_level_mp_term_name)
# summary(mod)
# 
# df_test <- aggregate(df2$percentage_change_abs, by=list(Category=df2$fun_Z_IMPC_percentile), FUN=mean)
# 
# plot(df_test$x ~ df_test$Category)
# mod <- lm(df_test$x ~ df_test$Category)
# summary(mod)
# 
# library(lme4)
# 
# library(lmerTest)
# 
# model = lmer(percentage_change_abs ~ fun_Z_IMPC_percentile + (1|top_level_mp_term_name),
#              data=df2,
#              REML=TRUE)
# anova(model)
# 
# MLexamp.6 <- lmer(percentage_change_abs ~ fun_Z_IMPC_percentile + (1|top_level_mp_term_name), data = df2)
# summary(MLexamp.6)
