### FIGURE 1 -- 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library("cowplot")

### FUNCTIONS

figure_1A <- function(M_df){

  M_plot <- ggplot(M_df, aes(x = n_total, y = n_functional)) +
  geom_point(aes(col=Classification)) +
  scale_color_manual(values=c('#CCCCCC', '#0000FF', '#FF0000')) +
  # ggtitle("Mouse") +
  xlab("Per-gene sum of all variant sites") +
  ylab('Per-gene sum of functional variant sites') +
  xlim(0, 300) +
  ylim(0, 150) +
  theme_bw() +
  theme(legend.position="none",
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1.5,1,1,1),"cm"))

return(M_plot)
}

figure_1B <- function(M_df){
  
  M_plot <- ggplot(M_df, aes(x = exp_fun, y = n_functional)) +
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
          plot.margin=unit(c(1.5,1,1,1),"cm"))
  return(M_plot)
}

# IMPORT
H_df <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")
M_df <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/M_funZ.csv")

# identify 5% extremes as most tolerant and intolerant
H_df$Classification <- "NA"
H_df$Classification[H_df$fun_Z_percentile_0.001 >= 0.98] <- "Most constrained"
H_df$Classification[H_df$fun_Z_percentile_0.001 <= 0.02] <- "Least constrained"
H_df$Classification <- as.factor(H_df$Classification)
H_df$Classification <- factor(H_df$Classification, levels = c('NA', "Least constrained", "Most constrained"))
H_df <- H_df[with(H_df, order(H_df$Classification)), ]

M_df$Classification <- "NA"
M_df$Classification[M_df$fun_Z_percentile_0.0001 >= 0.98] <- "Most constrained"
M_df$Classification[M_df$fun_Z_percentile_0.0001 <= 0.02] <- "Least constrained"
M_df$Classification <- as.factor(M_df$Classification)
M_df$Classification <- factor(M_df$Classification, levels = c('NA', "Least constrained", "Most constrained"))
M_df <- M_df[with(M_df, order(M_df$Classification)), ]

# remove extremes
# H_extremes <- subset(H_df, H_df$fun_Z_0.001 < -10 | H_df$fun_Z_0.001 > 10 )
# M_extremes <- subset(M_df, m_df$fun_Z_0.0001 < -10 | M_df$fun_Z_0.0001 > 10 )

### PLOT 
H_fig1A <- figure_1A(H_df)
M_fig1A <- figure_1A(M_df)
H_fig1B <- figure_1B(H_df)
M_fig1B <- figure_1B(M_df)

FigA <- plot_grid(H_fig1A, M_fig1A, ncol = 2, nrow = 1, labels = c("Human", "Mouse"), 
                 label_size = 20, hjust = -0.1, vjust = 2)
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_1A.jpg", FigA, ncol = 1, nrow = 1, base_height = 5, base_width = 10)

FigB <- plot_grid(M_fig1B, H_fig1B, ncol = 2, nrow = 1, labels = c("Mouse", "Human"), 
                  label_size = 20, hjust = -0.1, vjust = 2)
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_X_obs_exp.jpg", FigB, ncol = 1, nrow = 1, base_height = 5, base_width = 10)

Fig <- plot_grid(H_fig1A, M_fig1A, H_fig1B, M_fig1B, ncol = 2, nrow = 2, labels = c("Human", "Mouse", "Human", "Mouse"), 
                  label_size = 20, hjust = -0.1, vjust = 1.5)
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_1.jpg", Fig, ncol = 1, nrow = 1, base_height = 10, base_width = 10)



# x <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")
# y <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")
# 
# y <- y[,c("gene", "mis_z")]
# colnames(y) <- c("external_gene_name", "mis_z")
# z <- y[x, on = "external_gene_name"]
# z <- z[complete.cases(z),]
# plot(z$fun_Z_0.001, z$mis_z)
# cor.test(z$fun_Z_0.001, z$mis_z)
