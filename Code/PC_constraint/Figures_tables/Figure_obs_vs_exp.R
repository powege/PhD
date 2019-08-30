### FIGURE -- realtionship between observed and expected functional variants for MGP and 1KGP datasets

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library("cowplot")


# IMPORT
H_df <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF001_v2.csv")
M_df <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv")

# identify 5% extremes as most tolerant and intolerant

# calculate funZ percentile
percentile <- ecdf(H_df$fun_Z[!duplicated(H_df$external_gene_name)])
H_df$fun_Z_percentile <- percentile(H_df$fun_Z)
H_df$Classification <- "NA"
H_df$Classification[H_df$fun_Z_percentile >= 0.98] <- "Most constrained"
H_df$Classification[H_df$fun_Z_percentile <= 0.02] <- "Least constrained"
H_df$Classification <- as.factor(H_df$Classification)
H_df$Classification <- factor(H_df$Classification, levels = c('NA', "Least constrained", "Most constrained"))
H_df <- H_df[with(H_df, order(H_df$Classification)), ]

percentile <- ecdf(M_df$fun_Z[!duplicated(M_df$external_gene_name)])
M_df$fun_Z_percentile <- percentile(M_df$fun_Z)
M_df$Classification <- "NA"
M_df$Classification[M_df$fun_Z_percentile >= 0.98] <- "Most constrained"
M_df$Classification[M_df$fun_Z_percentile <= 0.02] <- "Least constrained"
M_df$Classification <- as.factor(M_df$Classification)
M_df$Classification <- factor(M_df$Classification, levels = c('NA', "Least constrained", "Most constrained"))
M_df <- M_df[with(M_df, order(M_df$Classification)), ]

# calculate totals 
# M_df$n_total <- M_df$n_missense + M_df$n_nonsense + M_df$n_synonymous
# H_df$n_total <- H_df$n_missense + H_df$n_nonsense + H_df$n_synonymous

# remove extremes
# H_extremes <- subset(H_df, H_df$fun_Z_0.001 < -10 | H_df$fun_Z_0.001 > 10 )
# M_extremes <- subset(M_df, m_df$fun_Z_0.0001 < -10 | M_df$fun_Z_0.0001 > 10 )

### PLOT 

M_plot <- ggplot(M_df, aes(x = exp_functional, y = n_functional)) +
  geom_point(aes(col=Classification)) +
  geom_abline(slope=1, intercept=0) +
  scale_color_manual(breaks = c("Least constrained", "Most constrained"),
                     values=c('gray80', 'blue', 'red')) +
  # ggtitle("Mouse") +
  xlab("Expected functional variant sites") +
  ylab('Observed functional variant sites') +
  xlim(0, 150) +
  ylim(0, 150) +
  theme_bw() +
  theme(
    legend.position="top",
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    text = element_text(size = 14),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_blank(),
    plot.margin=unit(c(1,1,1,1),"cm"))
M_plot

H_plot <- ggplot(H_df, aes(x = exp_functional, y = n_functional)) +
  geom_point(aes(col=Classification)) +
  geom_abline(slope=1, intercept=0) +
  scale_color_manual(breaks = c("Least constrained", "Most constrained"),
                     values=c('gray80', 'blue', 'red')) +
  # ggtitle("Mouse") +
  xlab("Expected functional variant sites") +
  ylab('Observed functional variant sites') +
  xlim(0, 150) +
  ylim(0, 150) +
  theme_bw() +
  theme(
    legend.position="top",
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    text = element_text(size = 14),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_blank(),
    plot.margin=unit(c(1,1,1,1),"cm"))
H_plot

H_total_syn <- sum(H_df$n_synonymous)
M_total_syn <- sum(M_df$n_synonymous)
H_total_mis <- sum(H_df$n_missense)
M_total_mis <- sum(M_df$n_missense)
H_total_non <- sum(H_df$n_nonsense)
M_total_non <- sum(M_df$n_nonsense)


Fig_out <- plot_grid(M_plot, H_plot, ncol = 2, nrow = 1, labels = c("Mouse", "Human"), 
                 label_size = 20, hjust = -0.1, vjust = 2)
save_plot("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_obs_exp.jpg", Fig_out, ncol = 1, nrow = 1, base_height = 5, base_width = 10)


# x <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")
# y <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")
# 
# y <- y[,c("gene", "mis_z")]
# colnames(y) <- c("external_gene_name", "mis_z")
# z <- y[x, on = "external_gene_name"]
# z <- z[complete.cases(z),]
# plot(z$fun_Z_0.001, z$mis_z)
# cor.test(z$fun_Z_0.001, z$mis_z)
