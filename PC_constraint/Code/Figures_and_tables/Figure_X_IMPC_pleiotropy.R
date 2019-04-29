### Relationship between constraint and number of MP terms (total, and top-level)
### numnber of MP terms adjusted for the number of phenotyping tests
### Prenatal phenotypes removed from analysis
### Data from IMPC release 9.1 (http://www.mousephenotype.org/data/release)

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(tidyr)
library("cowplot")

### FUNCTION for plotting lm equation on ggplot
linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 2),
            yy = format(abs(coef(k)[2]), digits = 2),
            r2 = format(summary(k)$r.squared, digits = 2));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)   
  }
  as.character(as.expression(eq));               
}

### FUNCTION for plotting r2 and p.val on ggplot
# k <- mod_all
r2p = function(k) {
  eq <- substitute(italic(r)^2~"="~rvalue*","~italic(p)~"="~pvalue, 
                   list(rvalue = format(summary(k)$r.squared, digits = 2), 
                        pvalue = format(anova(k)$'Pr(>F)'[1], digits = 2)))
    as.character(as.expression(eq))            
}

### IMPORT DATA
data <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")
impc <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/ALL_genotype_phenotype_release9.1.csv")
impc.all <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/IMPC_ALL_statistical_results.csv")

### FORMAT

### identify number of MP tests per KO 

no.mp.tests <- impc.all[,c("marker_symbol", "mp_term_name")]
no.mp.tests <- no.mp.tests[!is.na(no.mp.tests$mp_term_name),]
no.mp.tests <- no.mp.tests[no.mp.tests$mp_term_name != "",]
no.mp.tests <- separate_rows(no.mp.tests, mp_term_name, sep = ",")
no.mp.tests <- no.mp.tests[!duplicated(no.mp.tests),]
no.mp.tests <- as.data.table(table(no.mp.tests$marker_symbol))
colnames(no.mp.tests) <- c("external_gene_name", "n_MP_tests")

no.tlmp.tests <- impc.all[,c("marker_symbol", "top_level_mp_term_name")]
no.tlmp.tests <- no.tlmp.tests[!is.na(no.tlmp.tests$top_level_mp_term_name),]
no.tlmp.tests <- no.tlmp.tests[no.tlmp.tests$top_level_mp_term_name != "",]
no.tlmp.tests <- separate_rows(no.tlmp.tests, top_level_mp_term_name, sep = ",")
no.tlmp.tests <- no.tlmp.tests[!duplicated(no.tlmp.tests),]
no.tlmp.tests <- as.data.table(table(no.tlmp.tests$marker_symbol))
colnames(no.tlmp.tests) <- c("external_gene_name", "n_top_level_MP_tests")


### identify number of unique MPs per KO 

impc.mp <- impc[, c("marker_symbol", "mp_term_name")]
colnames(impc.mp) <- c("external_gene_name", "mp_term_name")
# remove genes with no top level phenotype term
impc.mp <- impc.mp[!is.na(impc.mp$mp_term_name),]
impc.mp <- impc.mp[impc.mp$mp_term_name != "",]
# split MPs into multiple categories
impc.mp <- separate_rows(impc.mp, mp_term_name, sep = ",")
# remove duplicates by gene name and mp term
impc.mp <- impc.mp[!duplicated(impc.mp[,c("external_gene_name", "mp_term_name")]),]
# count terms per gene
all.counts <- as.data.frame(table(impc.mp$external_gene_name))
colnames(all.counts) <- c("external_gene_name", "n_MP")

impc.tl <- impc[, c("marker_symbol", "top_level_mp_term_name")]
colnames(impc.tl) <- c("external_gene_name", "top_level_mp_term_name")
# remove genes with no top level phenotype term
impc.tl <- impc.tl[!is.na(impc.tl$top_level_mp_term_name),]
impc.tl <- impc.tl[impc.tl$top_level_mp_term_name != "",]
# split MPs into multiple categories
impc.tl <- separate_rows(impc.tl, top_level_mp_term_name, sep = ",")
# remove duplicates by gene name and mp term
impc.tl <- impc.tl[!duplicated(impc.tl[,c("external_gene_name", "top_level_mp_term_name")]),]
# count terms per gene
top.counts <- as.data.frame(table(impc.tl$external_gene_name))
colnames(top.counts) <- c("external_gene_name", "n_top_level_MP")


### idenmtify genes with no significant MP
sig_genes <- unique(impc$marker_symbol)
all_genes <- unique(impc.all$marker_symbol)
no_pheno_genes <- all_genes[-which(all_genes %in% sig_genes)]
no_pheno_genes <- data.frame(external_gene_name = no_pheno_genes,
                             n_MP = rep(0, length(no_pheno_genes)),
                             n_top_level_MP = rep(0, length(no_pheno_genes)))


# subset mouse funZ 
data <- data[,c("external_gene_name", "fun_Z_0.0001")]

# Merge data for plotting
df <- merge(no.mp.tests, no.tlmp.tests)
df <- merge(df, all.counts, all = T)
df <- merge(df, top.counts, all = T)
df <- merge(df, data)
# length(unique(df$external_gene_name))
df[is.na(df)] <- 0

# calculate funZ percentile for IMPC genes
percentile <- ecdf(df$fun_Z_0.0001)
df$fun_Z_IMPC_percentile <- percentile(df$fun_Z_0.0001)
df$fun_Z_IMPC_percentile <- ceiling((df$fun_Z_IMPC_percentile)*100)

# plot relationship between n_mp_tests and n_mps
plot(df$n_MP_tests, df$n_MP)
cor.test(df$n_MP_tests, df$n_MP)
plot(df$n_top_level_MP_tests, df$n_top_level_MP)
cor.test(df$n_top_level_MP_tests, df$n_top_level_MP)

# calculat hit rates
df$MP_hit_rate <- df$n_MP/df$n_MP_tests
df$top_level_MP_hit_rate <- df$n_top_level_MP/df$n_top_level_MP_tests

# count mean and median hit rate per percentile
df_mn_all <- aggregate(df$MP_hit_rate, by=list(Category=df$fun_Z_IMPC_percentile), FUN=mean)
df_mn_top <- aggregate(df$top_level_MP_hit_rate, by=list(Category=df$fun_Z_IMPC_percentile), FUN=mean)
df_md_all <- aggregate(df$MP_hit_rate, by=list(Category=df$fun_Z_IMPC_percentile), FUN=median)
df_md_top <- aggregate(df$top_level_MP_hit_rate, by=list(Category=df$fun_Z_IMPC_percentile), FUN=median)

# lm
mod_all <- lm(df_md_all$x~df_md_all$Category)
equation_all <- linear(mod_all)
r2p_all <- r2p(mod_all)
summary(mod_all)

plot_all <- ggplot(df_md_all, aes(x = Category, y = x)) +
  geom_point(alpha = 1/5) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "blue", size = 0.6, fullrange = T) +
  annotate("text", x = 37, y = max(df_md_all$x), label = r2p_all, colour="black", size = 5, parse=TRUE) +
  xlab("FunZ percentile bin") +
  ylab('Median MP term hit-rate') +
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
r2p_top <- r2p(mod_top)
summary(mod_top)

plot_top <- ggplot(df_md_top, aes(x = Category, y = x)) +
  geom_point(alpha = 1/5) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "blue", size = 0.6, fullrange = T) +
  annotate("text", x = 37, y = max(df_md_top$x), label = r2p_top, colour="black", size = 5, parse=TRUE) +
  xlab("FunZ percentile bin") +
  ylab('Median top-level MP term hit-rate') +
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
save_plot("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/Figure_X_pleiotropy.jpg", Fig, ncol = 1, nrow = 1, base_height = 5, base_width = 10)



#######

### EFFECT SIZE

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

####
# df <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/IMPC/IMPC_by_MP.csv")
# df <- subset(df, df$top_level_mp_term_name == "craniofacial phenotype (n=126)")
# cf <- unique(df$external_gene_name)
# cf <- impc[impc$marker_symbol %in% cf,]
# cf <- as.data.frame(table(cf$marker_symbol))
# colnames(cf) <- c("external_gene_name", "n_MP")
# cf <- df[cf, on = "external_gene_name"]


########

# gene list of pre-or-neonatal lethal knockouts
# lethal_genes <- impc[grep("lethal", impc$mp_term_name),]
# table(lethal_genes$mp_term_name)
# lethal_genes <- unique(lethal_genes$marker_symbol)
# embryonic_lethal <- unique(impc$marker_symbol[grep("embryonic lethal", impc$mp_term_name)])
# prenatal_lethal <- unique(impc$marker_symbol[grep("prenatal lethal", impc$mp_term_name)])
# prenatal_genes <- unique(c(embryonic_lethal, prenatal_lethal))

# rempve % from percentage_cjhange
# impc$percentage_change = substr(impc$percentage_change, 1, nchar(impc$percentage_change)-1)
# impc$percentage_change <- as.numeric(impc$percentage_change)
