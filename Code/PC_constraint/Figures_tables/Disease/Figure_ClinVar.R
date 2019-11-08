### SCRIPT that analyses ClinVar variant counts and gene constraint

rm(list = ls())
graphics.off()

library(data.table)
library("cowplot")

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

### FUNCTION for plotting Pearson's correlation on ggplot
spear_corr_eqn <- function(x, y, digits = 2) {
  corr_coef <- round(cor(x, y, method = "spearman"), digits = digits)
  paste("italic(rho) == ", corr_coef)
}


# IMPORT
CV <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_snps_QCed_VEP_v94.vcf", fill = T)
orth <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
H_con <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_MAF001_v2.csv")

### FORMAT

colnames(CV) <- c("CHROM", "POS", "ID", "REF", "ALT", "CAT", "FILTER", "Gene",
                  "Feature", "Feature_type", "Consequence", "IMPACT", "SYMBOL",
                  "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", 'CCDS')
CV <- CV[CV$Consequence %like% "synonymous_variant" |
           CV$Consequence %like% "stop_retained_variant" |
           CV$Consequence %like% "start_retained_variant" |
           CV$Consequence %like% "missense_variant" |
           CV$Consequence %like% "stop_gained" |
           CV$Consequence %like% "start_lost" |
           CV$Consequence %like% "stop_lost" |
           CV$Consequence %like% "splice_donor_variant" |
           CV$Consequence %like% "splice_acceptor_variant",]
test <- subset(CV, CV$Feature %in% orth$H_ensembl_transcript_id)
CV <- as.data.frame(table(CV$SYMBOL))
colnames(CV) <- c("H_external_gene_name", "n_ClinVar")

# number of variants
sum(CV$n_ClinVar)

# subset cds length
len <- H_con[,c("external_gene_name", "cds_length")]
colnames(len) <- paste0("H_", colnames(len))

# merge
df <- len[orth, on = "H_external_gene_name"]
df <- merge(df, CV, all = T)

# subset one 2 one ortholgs
# df <- df[df$orthology_type == "ortholog_one2one",]

# subset mouse genes
# m_df <- df[!duplicated(df$M_external_gene_name),]

# calculate one to one percentile
percentile <- ecdf(df$H_fun_Z[!duplicated(df$H_fun_Z)])
df$H_fun_Z_percentile <- percentile(df$H_fun_Z)
df$H_fun_Z_percentile <- ceiling((df$H_fun_Z_percentile)*100)
percentile <- ecdf(df$M_fun_Z[!duplicated(df$M_fun_Z)])
df$M_fun_Z_percentile <- percentile(df$M_fun_Z)
df$M_fun_Z_percentile <- ceiling((df$M_fun_Z_percentile)*100)

# calculate n_ClinVar per kb (n_ClinVar/cds_length)*1000
df$rel_ClinVar <- (df$n_ClinVar/df$H_cds_length)*1000
# length(unique(m_df$M_external_gene_name))

# calculate rel_ClinVar per gene
H_df <- df[,c("H_external_gene_name", "H_fun_Z_percentile", "rel_ClinVar")]
H_df$rel_ClinVar[is.na(H_df$rel_ClinVar)] <- 0
M_df <- df[,c("M_external_gene_name", "M_fun_Z_percentile", "rel_ClinVar")]
M_df$rel_ClinVar[is.na(M_df$rel_ClinVar)] <- 0

# number of genes
length(unique(H_df$H_external_gene_name))
length(unique(M_df$M_external_gene_name))

H_df_mn <- aggregate(H_df$rel_ClinVar, by=list(Category=H_df$H_fun_Z_percentile), FUN=mean)
H_df_md <- aggregate(H_df$rel_ClinVar, by=list(Category=H_df$H_fun_Z_percentile), FUN=median)
H_df_sum <- aggregate(H_df$rel_ClinVar, by=list(Category=H_df$H_fun_Z_percentile), FUN=sum)

M_df_mn <- aggregate(M_df$rel_ClinVar, by=list(Category=M_df$M_fun_Z_percentile), FUN=mean)
M_df_md <- aggregate(M_df$rel_ClinVar, by=list(Category=M_df$M_fun_Z_percentile), FUN=median)
M_df_sum <- aggregate(M_df$rel_ClinVar, by=list(Category=M_df$M_fun_Z_percentile), FUN=sum)


# fit model
cor.test(H_df_mn$x, H_df_mn$Category, method = "spearman")
H_mn_mod <- lm(x ~ Category, data = H_df_mn)
summary(H_mn_mod)
H_lm_equation <- linear(H_mn_mod)

cor.test(M_df_mn$x, M_df_mn$Category, method = "spearman")
M_mn_mod <- lm(x ~ Category, data = M_df_mn)
summary(M_mn_mod)

### PLOT 

H_lm_equation <- linear(H_mn_mod)
M_lm_equation <- linear(M_mn_mod)
H_spearman <- spear_corr_eqn(H_df_mn$x, H_df_mn$Category, digits = 3)
M_spearman <- spear_corr_eqn(M_df_mn$x, M_df_mn$Category, digits = 3)

M_plot <- ggplot() +
  geom_bar(data=M_df_md, aes(y=x,x=Category), color = "red", fill="red", stat="identity", width = 0.75) + 
  geom_point(data=M_df_mn, aes(y=x, x=Category), alpha = 1/5) +
  geom_smooth(data=M_df_mn, aes(y=x, x=Category), color = "black", method='lm', size = 0.6, fullrange = T, linetype = "longdash", se = F) +
  # scale_color_gradient(low="blue", high="red") +
  annotate("text", x = 20, y = 3.5, label = M_spearman, colour="black", size = 5, parse=TRUE) +
  xlab("Mouse funZ percentile bin") +
  ylab('Pathogenic variants per kb\nin human orthologue') +ylim(0, 4) +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1.5,1,1,1),"cm"))
M_plot

H_plot <- ggplot() +
  geom_bar(data=H_df_md, aes(y=x,x=Category), color = "red", fill="red", stat="identity", width = 0.75) + 
  geom_point(data=H_df_mn, aes(y=x, x=Category), alpha = 1/5) +
  geom_smooth(data=H_df_mn, aes(y=x, x=Category), color = "black", method='lm', size = 0.6, fullrange = T, linetype = "longdash", se = F) +
  # scale_color_gradient(low="blue", high="red") +
  annotate("text", x = 20, y = 3.5, label = H_spearman, colour="black", size = 5, parse=TRUE) +
  xlab("Human funZ percentile bin") +
  ylab('Pathogenic variants per kb') +ylim(0, 4) +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1.5,1,1,1),"cm"))
H_plot

Fig <- plot_grid(M_plot, H_plot, ncol = 2, nrow = 1, labels = c("A", "B"), 
                 label_size = 20, hjust = -0.1, vjust = 1.5)
Fig
save_plot("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_ClinVar.jpg", Fig, ncol = 1, nrow = 1, base_height = 5, base_width = 10)


#####

# ### PLOT WITH DIFFERENT CONSTRAINT METRIC
# 
# # IMPORT
# CV <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_pathogenic_snps_QCed_VEP_v94.vcf", fill = T)
# orth <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
# H_con <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_MAF001_v2.csv")
# 
# ### FORMAT
# 
# colnames(CV) <- c("CHROM", "POS", "ID", "REF", "ALT", "CAT", "FILTER", "Gene",
#                   "Feature", "Feature_type", "Consequence", "IMPACT", "SYMBOL",
#                   "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", 'CCDS')
# CV <- CV[CV$Consequence %like% "synonymous_variant" |
#            CV$Consequence %like% "stop_retained_variant" |
#            CV$Consequence %like% "start_retained_variant" |
#            CV$Consequence %like% "missense_variant" |
#            CV$Consequence %like% "stop_gained" |
#            CV$Consequence %like% "start_lost" |
#            CV$Consequence %like% "stop_lost" |
#            CV$Consequence %like% "splice_donor_variant" |
#            CV$Consequence %like% "splice_acceptor_variant",]
# test <- subset(CV, CV$Feature %in% orth$H_ensembl_transcript_id)
# CV <- as.data.frame(table(CV$SYMBOL))
# colnames(CV) <- c("H_external_gene_name", "n_ClinVar")
# 
# # number of variants
# sum(CV$n_ClinVar)
# 
# # subset cds length
# len <- H_con[,c("external_gene_name", "cds_length")]
# colnames(len) <- paste0("H_", colnames(len))
# 
# # merge
# df <- len[orth, on = "H_external_gene_name"]
# df <- merge(df, CV, all = T)
# 
# # subset one 2 one ortholgs
# # df <- df[df$orthology_type == "ortholog_one2one",]
# 
# # subset mouse genes
# # m_df <- df[!duplicated(df$M_external_gene_name),]
# 
# # df$H_fun_Z <- df$H_OE_ratio
# # df$M_fun_Z <- df$M_OE_ratio
# df$H_fun_Z <- df$H_fun_Z_adj
# df$M_fun_Z <- df$M_fun_Z_adj
# 
# # calculate one to one percentile
# percentile <- ecdf(df$H_fun_Z[!duplicated(df$H_fun_Z)])
# df$H_fun_Z_percentile <- percentile(df$H_fun_Z)
# df$H_fun_Z_percentile <- ceiling((df$H_fun_Z_percentile)*100)
# percentile <- ecdf(df$M_fun_Z[!duplicated(df$M_fun_Z)])
# df$M_fun_Z_percentile <- percentile(df$M_fun_Z)
# df$M_fun_Z_percentile <- ceiling((df$M_fun_Z_percentile)*100)
# 
# # calculate n_ClinVar per kb (n_ClinVar/cds_length)*1000
# df$rel_ClinVar <- (df$n_ClinVar/df$H_cds_length)*1000
# # length(unique(m_df$M_external_gene_name))
# 
# # calculate rel_ClinVar per gene
# H_df <- df[,c("H_external_gene_name", "H_fun_Z_percentile", "rel_ClinVar")]
# H_df$rel_ClinVar[is.na(H_df$rel_ClinVar)] <- 0
# M_df <- df[,c("M_external_gene_name", "M_fun_Z_percentile", "rel_ClinVar")]
# M_df$rel_ClinVar[is.na(M_df$rel_ClinVar)] <- 0
# 
# # number of genes
# length(unique(H_df$H_external_gene_name))
# length(unique(M_df$M_external_gene_name))
# 
# H_df_mn <- aggregate(H_df$rel_ClinVar, by=list(Category=H_df$H_fun_Z_percentile), FUN=mean)
# H_df_md <- aggregate(H_df$rel_ClinVar, by=list(Category=H_df$H_fun_Z_percentile), FUN=median)
# H_df_sum <- aggregate(H_df$rel_ClinVar, by=list(Category=H_df$H_fun_Z_percentile), FUN=sum)
# 
# M_df_mn <- aggregate(M_df$rel_ClinVar, by=list(Category=M_df$M_fun_Z_percentile), FUN=mean)
# M_df_md <- aggregate(M_df$rel_ClinVar, by=list(Category=M_df$M_fun_Z_percentile), FUN=median)
# M_df_sum <- aggregate(M_df$rel_ClinVar, by=list(Category=M_df$M_fun_Z_percentile), FUN=sum)
# 
# 
# # fit model
# cor.test(H_df_mn$x, H_df_mn$Category, method = "spearman")
# H_mn_mod <- lm(x ~ Category, data = H_df_mn)
# summary(H_mn_mod)
# H_lm_equation <- linear(H_mn_mod)
# 
# cor.test(M_df_mn$x, M_df_mn$Category, method = "spearman")
# M_mn_mod <- lm(x ~ Category, data = M_df_mn)
# summary(M_mn_mod)
# 
# ### PLOT
# 
# H_lm_equation <- linear(H_mn_mod)
# M_lm_equation <- linear(M_mn_mod)
# H_spearman <- spear_corr_eqn(H_df_mn$x, H_df_mn$Category, digits = 3)
# M_spearman <- spear_corr_eqn(M_df_mn$x, M_df_mn$Category, digits = 3)
# 
# M_plot <- ggplot() +
#   geom_bar(data=M_df_md, aes(y=x,x=Category), color = "red", fill="red", stat="identity", width = 0.75) +
#   geom_point(data=M_df_mn, aes(y=x, x=Category), alpha = 1/5) +
#   geom_smooth(data=M_df_mn, aes(y=x, x=Category), color = "blue", method='lm') +
#   # scale_color_gradient(low="blue", high="red") +
#   annotate("text", x = 20, y = 3.5, label = M_spearman, colour="black", size = 5, parse=TRUE) +
#   xlab("Mouse funZ percentile bin") +
#   ylab('Pathogenic variants per kb\nin human orthologue') +ylim(0, 4) +
#   theme_bw() +
#   theme(legend.position="none",
#         text = element_text(size = 14),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.background=element_blank(),
#         plot.margin=unit(c(1.5,1,1,1),"cm"))
# M_plot
# 
# H_plot <- ggplot() +
#   geom_bar(data=H_df_md, aes(y=x,x=Category), color = "red", fill="red", stat="identity", width = 0.75) +
#   geom_point(data=H_df_mn, aes(y=x, x=Category), alpha = 1/5) +
#   geom_smooth(data=H_df_mn, aes(y=x, x=Category), color = "blue", method='lm') +
#   # scale_color_gradient(low="blue", high="red") +
#   annotate("text", x = 20, y = 3.5, label = H_spearman, colour="black", size = 5, parse=TRUE) +
#   xlab("Human funZ percentile bin") +
#   ylab('Pathogenic variants per kb') +ylim(0, 4) +
#   theme_bw() +
#   theme(legend.position="none",
#         text = element_text(size = 14),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.background=element_blank(),
#         plot.margin=unit(c(1.5,1,1,1),"cm"))
# H_plot
# 
# Fig <- plot_grid(M_plot, H_plot, ncol = 2, nrow = 1, labels = c("A", "B"),
#                  label_size = 20, hjust = -0.1, vjust = 1.5)
# Fig

# # Calculate fold change relative to 1st percentile
# colnames(df_mn) <- c("percentile", "mean_CV_per_kb")
# df_mn$fold_change <- df_mn$mean_CV_per_kb / df_mn$mean_CV_per_kb[df_mn$percentile == 1]
# 
# plot(df_mn$mean_CV_per_kb ~ df_mn$percentile)
# plot(df_mn$fold_change ~ df_mn$percentile)
# 
# gobbler <- function(h_gerp){
#   percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                         "51-75", "76-90", "91-95", "96-98", "99", "100")
# 
#   mean_fold_change <- rep(NA, 12)
#   mean_fold_change[1] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(1)])
#   mean_fold_change[2] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(2)])
#   mean_fold_change[3] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])
#   mean_fold_change[4] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])
#   mean_fold_change[5] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])
#   mean_fold_change[6] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])
#   mean_fold_change[7] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])
#   mean_fold_change[8] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])
#   mean_fold_change[9] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])
#   mean_fold_change[10] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])
#   mean_fold_change[11] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(99)])
#   mean_fold_change[12] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(100)])
# 
#   upper_range_fold_change <- rep(NA, 12)
#   upper_range_fold_change[1] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(1)])[1]
#   upper_range_fold_change[2] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(2)])[1]
#   upper_range_fold_change[3] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])[1]
#   upper_range_fold_change[4] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])[1]
#   upper_range_fold_change[5] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])[1]
#   upper_range_fold_change[6] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])[1]
#   upper_range_fold_change[7] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])[1]
#   upper_range_fold_change[8] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])[1]
#   upper_range_fold_change[9] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])[1]
#   upper_range_fold_change[10] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])[1]
#   upper_range_fold_change[11] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(99)])[1]
#   upper_range_fold_change[12] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(100)])[1]
# 
#   lower_range_fold_change <- rep(NA, 12)
#   lower_range_fold_change[1] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(1)])[2]
#   lower_range_fold_change[2] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(2)])[2]
#   lower_range_fold_change[3] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])[2]
#   lower_range_fold_change[4] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])[2]
#   lower_range_fold_change[5] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])[2]
#   lower_range_fold_change[6] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])[2]
#   lower_range_fold_change[7] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])[2]
#   lower_range_fold_change[8] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])[2]
#   lower_range_fold_change[9] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])[2]
#   lower_range_fold_change[10] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])[2]
#   lower_range_fold_change[11] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(99)])[2]
#   lower_range_fold_change[12] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(100)])[2]
# 
#   out <- data.frame(percentile_group = percentile_group,
#                     rank = c(1:12),
#                     mean_fold_change = mean_fold_change,
#                     lower_range_fold_change = lower_range_fold_change,
#                     upper_range_fold_change = upper_range_fold_change
#                     # species = rep(h_gerp$species[1], 12),
#                     # annotation = rep(h_gerp$annotation[1], 12)
#                     )
#   return(out)
# }
# 
# gobbler2 <- function(h_gerp){
#   percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                         "51-75", "76-90", "91-95", "96-98", "99", "100")
# 
#   mean_fold_change <- rep(NA, 12)
#   mean_fold_change[1] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(1)])
#   mean_fold_change[2] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(2)])
#   mean_fold_change[3] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])
#   mean_fold_change[4] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])
#   mean_fold_change[5] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])
#   mean_fold_change[6] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])
#   mean_fold_change[7] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])
#   mean_fold_change[8] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])
#   mean_fold_change[9] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])
#   mean_fold_change[10] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])
#   mean_fold_change[11] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(99)])
#   mean_fold_change[12] <- median(h_gerp$fold_change[h_gerp$percentile %in% c(100)])
# 
#   upper_quantile_fold_change <- rep(NA, 12)
#   upper_quantile_fold_change[1] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(1)], 0.75)
#   upper_quantile_fold_change[2] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(2)], 0.75)
#   upper_quantile_fold_change[3] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)], 0.75)
#   upper_quantile_fold_change[4] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)], 0.75)
#   upper_quantile_fold_change[5] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)], 0.75)
#   upper_quantile_fold_change[6] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)], 0.75)
#   upper_quantile_fold_change[7] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)], 0.75)
#   upper_quantile_fold_change[8] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)], 0.75)
#   upper_quantile_fold_change[9] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)], 0.75)
#   upper_quantile_fold_change[10] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)], 0.75)
#   upper_quantile_fold_change[11] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(99)], 0.75)
#   upper_quantile_fold_change[12] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(100)], 0.75)
# 
#   lower_quantile_fold_change <- rep(NA, 12)
#   lower_quantile_fold_change[1] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(1)], 0.25)
#   lower_quantile_fold_change[2] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(2)], 0.25)
#   lower_quantile_fold_change[3] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)], 0.25)
#   lower_quantile_fold_change[4] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)], 0.25)
#   lower_quantile_fold_change[5] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)], 0.25)
#   lower_quantile_fold_change[6] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)], 0.25)
#   lower_quantile_fold_change[7] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)], 0.25)
#   lower_quantile_fold_change[8] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)], 0.25)
#   lower_quantile_fold_change[9] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)], 0.25)
#   lower_quantile_fold_change[10] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)], 0.25)
#   lower_quantile_fold_change[11] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(99)], 0.25)
#   lower_quantile_fold_change[12] <- quantile(h_gerp$fold_change[h_gerp$percentile %in% c(100)], 0.25)
# 
#   out <- data.frame(percentile_group = percentile_group,
#                     rank = c(1:12),
#                     mean_fold_change = mean_fold_change,
#                     lower_quantile_fold_change = lower_quantile_fold_change,
#                     upper_quantile_fold_change = upper_quantile_fold_change
#                     # species = rep(h_gerp$species[1], 12),
#                     # annotation = rep(h_gerp$annotation[1], 12)
#   )
#   return(out)
# }
# 
# test <- gobbler2(df_mn)
# p_h <- ggplot(test, aes(x=rank, y=mean_fold_change)) +
#   # geom_point() +
#   geom_line(size=1.5) +
#   geom_errorbar(aes(ymax = upper_quantile_fold_change, ymin = lower_quantile_fold_change),
#                 stat = "identity",
#                 width=0.2,
#                 size=1) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Pathogenic variants per kb\n(fold change versus 1st percentile)") +
#   ggtitle("Human") +
#   # scale_y_continuous( position = "right",
#   #                     trans = log2_trans(),
#   #                     # breaks = trans_breaks("log2", function(x) 2^x)) +
#   #                     breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40),
#   #                     limits = c(0.05, 40)) +
#   scale_x_continuous(breaks=c(1:12),
#                      # trans = "reverse",
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     # axis.title.y = element_blank(),
#     legend.title = element_blank(),
#     legend.justification=c(1,0),
#     legend.position=c(0.23, 0.70),
#     legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# p_h
# 
# 
# Fig <- plot_grid(m_plot, h_plot, ncol = 2, nrow = 1, labels = c("A", "B"),
#                  label_size = 20, hjust = -0.1, vjust = 2)
# Fig

