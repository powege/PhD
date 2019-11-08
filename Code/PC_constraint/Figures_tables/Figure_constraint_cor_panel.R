rm(list=ls())
graphics.off()

library(corrplot)
library(tidyr)
library(RColorBrewer)
library("cowplot")
library(data.table)

### FUNCTION for plotting Spearman's correlation on ggplot
rs_corr_eqn <- function(x,y, digits) {
  corr_coef <- round(cor(x, y, method = "spearman"), digits = digits)
  paste("italic(rho) == ", corr_coef)
}

### FUNCTION for plotting Spearman's correlation on ggplot
r_corr_eqn <- function(x,y, digits) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}

### IMPORT

rvis <- fread("~/Dropbox/PhD/Data/Constraint/RVIS_Petrovski2013_ESP.csv")
gnomad <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt")
funZ <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")

### FORMAT

gnomad <- gnomad[,c("gene", "mis_z", "pLI", "oe_lof_upper")]
colnames(gnomad) <- c("H_external_gene_name", "mis_z", "pLI", "oe_lof_upper")
rvis <- rvis[,c("HGNC gene", "Residual Variation Intolerance Score")]
colnames(rvis) <- c("H_external_gene_name", "RVIS")
funZ <- funZ[,c("H_external_gene_name", "M_external_gene_name", 
                "M_fun_Z", "H_fun_Z", 
                "M_OE_ratio", "H_OE_ratio", 
                "M_fun_Z_MSE", "H_fun_Z_MSE", 
                # "M_fun_Z_adj", "H_fun_Z_adj", 
                "orthology_type")]

# subset one to one orthologues
funZ_o2o <- subset(funZ, funZ$orthology_type == "ortholog_one2one")

# merge
dt_all <- gnomad[rvis, on = "H_external_gene_name"]
dt_all <- dt_all[funZ_o2o, on = "H_external_gene_name"]
dt_all <- dt_all[!duplicated(dt_all),]
dt_all <- dt_all[complete.cases(dt_all),]

# calculate percentiles
# percentile <- ecdf(dt_all$mis_z[!duplicated(dt_all$H_external_gene_name)])
# dt_all$mis_z_percentile <- percentile(dt_all$mis_z)
# percentile <- ecdf(dt_all$pLI[!duplicated(dt_all$H_external_gene_name)])
# dt_all$pLI_percentile <- percentile(dt_all$pLI)
# percentile <- ecdf(dt_all$oe_lof_upper[!duplicated(dt_all$H_external_gene_name)])
# dt_all$oe_lof_upper_percentile <- percentile(dt_all$oe_lof_upper)
# percentile <- ecdf(dt_all$RVIS[!duplicated(dt_all$H_external_gene_name)])
# dt_all$RVIS_percentile <- percentile(dt_all$RVIS)
# percentile <- ecdf(dt_all$M_fun_Z[!duplicated(dt_all$M_external_gene_name)])
# dt_all$M_fun_Z_percentile <- percentile(dt_all$M_fun_Z)
# percentile <- ecdf(dt_all$H_fun_Z[!duplicated(dt_all$H_external_gene_name)])
# dt_all$H_fun_Z_percentile <- percentile(dt_all$H_fun_Z)
# percentile <- ecdf(dt_all$M_OE_ratio[!duplicated(dt_all$M_external_gene_name)])
# dt_all$M_OE_ratio_percentile <- percentile(dt_all$M_OE_ratio)
# percentile <- ecdf(dt_all$H_OE_ratio[!duplicated(dt_all$H_external_gene_name)])
# dt_all$H_OE_ratio_percentile <- percentile(dt_all$H_OE_ratio)
# percentile <- ecdf(dt_all$M_fun_Z_adj[!duplicated(dt_all$M_external_gene_name)])
# dt_all$M_fun_Z_adj_percentile <- percentile(dt_all$M_fun_Z_adj)
# percentile <- ecdf(dt_all$H_fun_Z_adj[!duplicated(dt_all$H_external_gene_name)])
# dt_all$H_fun_Z_adj_percentile <- percentile(dt_all$H_fun_Z_adj)

# reverse RVIS and oe_lof_upper
# dt_all$RVIS_percentile <- 0-dt_all$RVIS_percentile
# dt_all$oe_lof_upper_percentile <- 0-dt_all$oe_lof_upper_percentile

# correlation matrix
# dt_plot <- cor(dt_all[,c("M_fun_Z_percentile",
#                          "H_fun_Z_percentile",
#                          # "M_fun_Z_adj_percentile",
#                          # "H_fun_Z_adj_percentile",
#                          # "M_OE_ratio_percentile",
#                          # "H_OE_ratio_percentile",
#                          "RVIS_percentile",
#                          "mis_z_percentile",
#                          "pLI_percentile",
#                          "oe_lof_upper_percentile")], method = "spearman")
dt_plot <- cor(dt_all[,c("M_fun_Z",
                         "H_fun_Z",
                         "RVIS",
                         "mis_z",
                         "pLI",
                         "oe_lof_upper")], method = "spearman")
dt_plot <- as.matrix(dt_plot)
colnames(dt_plot) <- c("Mouse funZ", 
                       "Human funZ",
                       # "Mouse funZ (adjusted)", 
                       # "Human funZ (adjusted)", 
                       # "Mouse OE ratio",
                       # "Human OE ratio",
                       "RVIS", 
                       "missense z-score", 
                       "pLI", 
                       "LOEUF")
rownames(dt_plot) <- c("Mouse funZ", 
                       "Human funZ",
                       # "Mouse funZ (adjusted)", 
                       # "Human funZ (adjusted)", 
                       # "Mouse OE ratio",
                       # "Human OE ratio",
                       "RVIS", 
                       "missense z-score", 
                       "pLI", 
                       "LOEUF")
melted_cormat <- melt(dt_plot, na.rm = TRUE)
melted_cormat$value <- round(melted_cormat$value, digits = 2)

dt_plot1 <- funZ[,c("M_external_gene_name", "H_external_gene_name", "M_fun_Z", "H_fun_Z", "orthology_type")]

# caclulate percentile for orhtologues
percentile <- ecdf(dt_plot1$M_fun_Z[!duplicated(dt_plot1$M_external_gene_name)]) # remove duplicates 
dt_plot1$M_fun_Z_percentile <- percentile(dt_plot1$M_fun_Z)
percentile <- ecdf(dt_plot1$H_fun_Z[!duplicated(dt_plot1$H_external_gene_name)]) # remove duplicates
dt_plot1$H_fun_Z_percentile <- percentile(dt_plot1$H_fun_Z)

### identify overlapping extremes as most tolerant and intolerant
dt_plot1$HM_funZ_class <- "NA"
dt_plot1$HM_funZ_class[dt_plot1$H_fun_Z_percentile >= 0.9 & dt_plot1$M_fun_Z_percentile >= 0.9] <- "10% most constrained in both species"

dt_plot1$HM_funZ_class <- as.factor(dt_plot1$HM_funZ_class)
dt_plot1$HM_funZ_class <- factor(dt_plot1$HM_funZ_class, levels = c("NA", "10% most constrained in both species"))
table(dt_plot1$HM_funZ_class)

### COR TEST
cor.test(dt_plot1$M_fun_Z, dt_plot1$H_fun_Z)

### PLOT 1 

cor_text <- r_corr_eqn(dt_plot1$H_fun_Z, dt_plot1$M_fun_Z, digits = 2)
p1 <- ggplot(dt_plot1, aes(x = M_fun_Z, y = H_fun_Z)) +
  geom_point(aes(col=HM_funZ_class)) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T, linetype = "longdash") +
  annotate("text", x = -7, y = 8, label = cor_text, colour="black", size = 5, parse=TRUE) +
  scale_color_manual(breaks = c("10% most constrained in both species"), values=c('gray80', 'red')) +
  xlab("Mouse funZ") +
  ylab('Human funZ') +
  xlim(-10, 10) +
  ylim(-10, 10) +
  theme_bw() +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
p1

### PLOT 2

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red3", high = "blue3", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Rank\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(vjust = 1, 
                                   size = 14, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value, fontface=2), color = "black", size = 5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(1.3, 0.0),
    # legend.direction = "vertical"
    text = element_text(size = 14)
    )+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 12,
                               title.position = "top", title.hjust = 0.5))
ggheatmap 

### COMBINE 

Fig_out <- plot_grid(p1, ggheatmap, ncol = 2, nrow = 1, labels = c("A", "B"), 
                     label_size = 20, hjust = -0.1, vjust = 2)
Fig_out
save_plot("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_HM_cor.jpg", Fig_out, ncol = 1, nrow = 1, base_height = 6, base_width = 12)


