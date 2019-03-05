rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

df <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")

table(df$orthology_type)
o2o <- subset(df, df$orthology_type == "ortholog_one2one") # one mouse to one human
m2m <- subset(df, df$orthology_type == "ortholog_many2many") # many mouse to many human
o2m <- subset(df, df$orthology_type == "ortholog_one2many") # one human to many mouse

# MOUSE df
M_df <- df[,c("M_external_gene_name", "orthology_type", "M_fun_Z_0.0001", "M_fun_Z_percentile_0.0001")]
M_df <- M_df[!duplicated(M_df$M_external_gene_name),]
M_df <- M_df[complete.cases(M_df),]

# calculate funZ percentile for IMPC genes
percentile <- ecdf(M_df$M_fun_Z_0.0001)
M_df$M_fun_Z_percentile <- percentile(M_df$M_fun_Z_0.0001)
# plot(M_df$M_fun_Z_percentile_0.0001, M_df$M_fun_Z_percentile)

# N 
# M_o2o <- M_df[M_df$orthology_type == "ortholog_one2one",]
# M_o2m <- M_df[M_df$orthology_type == "ortholog_one2many",]
# M_m2m <- M_df[M_df$orthology_type == "ortholog_many2many",]
M_n.o2o <- length(unique(M_df$M_external_gene_name[M_df$orthology_type == "ortholog_one2one"]))
M_n.o2m <- length(unique(M_df$M_external_gene_name[M_df$orthology_type == "ortholog_one2many"]))
M_n.m2m <- length(unique(M_df$M_external_gene_name[M_df$orthology_type == "ortholog_many2many"]))


M_df$orthology_type[M_df$orthology_type == "ortholog_one2one"] <- paste0("one-to-one\n(n=", M_n.o2o, ")")
M_df$orthology_type[M_df$orthology_type == "ortholog_one2many"] <- paste0("one-to-many\n(n=", M_n.o2m, ")")
M_df$orthology_type[M_df$orthology_type == "ortholog_many2many"] <- paste0("many-to-many\n(n=", M_n.m2m, ")")

M_df$orthology_type <- reorder(M_df$orthology_type, M_df$M_fun_Z_percentile, FUN = median, order = T)
table(M_df$orthology_type)
M_comparisons <- list( c("one-to-one\n(n=13787)", "one-to-many\n(n=1302)"), 
                        c("one-to-one\n(n=13787)", "many-to-many\n(n=333)"),
                        c("many-to-many\n(n=333)", "one-to-many\n(n=1302)"))

M_box <- ggplot(M_df, aes(x=orthology_type, y=M_fun_Z_percentile, fill = orthology_type)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = M_comparisons,
    method = "wilcox.test", 
    label = "p.signif",
    label.y = c(1.15, 1.25, 1.05),
    paired = F
    # size = 3,
  ) +
  ylab("Mouse funZ percentile rank") +
  xlab("") +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1.5,1,1,1),"cm")) 
M_box


# HUMAN df
H_df <- df[,c("H_external_gene_name", "orthology_type", "H_fun_Z_0.001", "H_fun_Z_percentile_0.001")]
H_df <- H_df[!duplicated(H_df$H_external_gene_name),]
H_df <- H_df[complete.cases(H_df),]

# calculate funZ percentile for IMPC genes
percentile <- ecdf(H_df$H_fun_Z_0.001)
H_df$H_fun_Z_percentile <- percentile(H_df$H_fun_Z_0.001)
# plot(H_df$H_fun_Z_percentile_0.001, H_df$H_fun_Z_percentile)

# N 
# H_o2o <- H_df[H_df$orthology_type == "ortholog_one2one",]
# H_o2m <- H_df[H_df$orthology_type == "ortholog_one2many",]
# H_m2m <- H_df[H_df$orthology_type == "ortholog_many2many",]
H_n.o2o <- length(unique(H_df$H_external_gene_name[H_df$orthology_type == "ortholog_one2one"]))
H_n.o2m <- length(unique(H_df$H_external_gene_name[H_df$orthology_type == "ortholog_one2many"]))
H_n.m2m <- length(unique(H_df$H_external_gene_name[H_df$orthology_type == "ortholog_many2many"]))


H_df$orthology_type[H_df$orthology_type == "ortholog_one2one"] <- paste0("one-to-one\n(n=", H_n.o2o, ")")
H_df$orthology_type[H_df$orthology_type == "ortholog_one2many"] <- paste0("one-to-many\n(n=", H_n.o2m, ")")
H_df$orthology_type[H_df$orthology_type == "ortholog_many2many"] <- paste0("many-to-many\n(n=", H_n.m2m, ")")

H_df$orthology_type <- reorder(H_df$orthology_type, H_df$H_fun_Z_percentile, FUN = median, order = T)
table(H_df$orthology_type)
H_comparisons <- list( c("one-to-one\n(n=13787)", "one-to-many\n(n=902)"), 
                       c("one-to-one\n(n=13787)", "many-to-many\n(n=293)"),
                       c("many-to-many\n(n=293)", "one-to-many\n(n=902)"))

H_box <- ggplot(H_df, aes(x=orthology_type, y=H_fun_Z_percentile, fill = orthology_type)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = H_comparisons,
    method = "wilcox.test", 
    label = "p.signif",
    label.y = c(1.15, 1.25, 1.05),
    paired = F
    # size = 3,
  ) +
  ylab("Human funZ percentile rank") +
  xlab("") +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1.5,1,1,1),"cm")) 
H_box

Fig <- plot_grid(H_box, M_box, ncol = 2, nrow = 1, labels = c("Human", "Mouse"), 
                  label_size = 20, hjust = -0.1, vjust = 2)
Fig
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_X_ortholog_box.jpg", Fig, ncol = 1, nrow = 1, base_height = 5.5, base_width = 11)



