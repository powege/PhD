# http://www.informatics.jax.org/disease/DOID:4

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(ggpubr)
library("cowplot")


### SET RELATIVE PATHS
mgi.path <- "~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/"

### IMPORT DATA
newC <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
RVIS_esp <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/RVIS_Petrovski2013_ESP.csv")
Lek <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")
mgi.AE <- fread(paste0(mgi.path, "MGIDO_anatomical_entity.txt"))
mgi.CP <- fread(paste0(mgi.path, "MGIDO_cellular_proliferation.txt"))
mgi.GD <- fread(paste0(mgi.path, "MGIDO_genetic_disease.txt"))
mgi.ID <- fread(paste0(mgi.path, "MGIDO_infectious_agent.txt"))
mgi.MH <- fread(paste0(mgi.path, "MGIDO_mental_health.txt"))
mgi.M <- fread(paste0(mgi.path, "MGIDO_metabolism.txt"))
mgi.PD <- fread(paste0(mgi.path, "MGIDO_physical_disorder.txt"))
mgi.S <- fread(paste0(mgi.path, "MGIDO_syndrome.txt"))


### FORMAT

# colnames
colnames(RVIS_esp) <- c("H_external_gene_name", "RVIS_ESP", "RVIS_ESP_percentile")
Lek <- Lek[,c("gene", "mis_z", "lof_z", "pLI")]
colnames(Lek) <- c("H_external_gene_name", "mis_z_Lek", "lof_z_Lek", "pLI_Lek")

# percentiles
RVIS_esp$RVIS_ESP_percentile <- RVIS_esp$RVIS_ESP_percentile/100
percentile <- ecdf(Lek$mis_z_Lek)
Lek$mis_z_percentile_Lek <- percentile(Lek$mis_z_Lek)
percentile <- ecdf(Lek$lof_z_Lek)
Lek$lof_z_percentile_Lek <- percentile(Lek$lof_z_Lek)
percentile <- ecdf(Lek$pLI_Lek)
Lek$pLI_percentile_Lek <- percentile(Lek$pLI_Lek)

df <- merge(Lek, RVIS_esp, all = T)
df <- merge(df, newC, all = T)
df <- df[!is.na(df$H_external_gene_name),]
df <- df[df$H_external_gene_name != "",]

# Set gene lists

# add disease catagory
mgi.AE$Disease_cat <- "anatomical_entity"
mgi.CP$Disease_cat <- "cellular_proliferation"
mgi.GD$Disease_cat <- "genetic_disease"
mgi.M$Disease_cat <- "metabolism"
mgi.MH$Disease_cat <- "mental_health"
mgi.PD$Disease_cat <- "physical_disorder"
mgi.S$Disease_cat <- "syndrome"
mgi.ID$Disease_cat <- "infectious_agent"

# combine MGI
mgi <- rbind(mgi.AE, mgi.CP, mgi.GD, mgi.ID, mgi.M, mgi.MH, mgi.PD, mgi.S)

colnames(mgi) <- c("Gene_Catagory", "Disease_term", "M_external_gene_name", "H_external_gene_name", 
                   "Other_features", "Mouse_models", "Homology_source", "BLANK", "Disease_cat")
mgi <- mgi[,c("Gene_Catagory", "Disease_term", "M_external_gene_name", "H_external_gene_name", "Disease_cat")]

# remove all *
mgi$M_external_gene_name <- gsub("[*]", "", mgi$M_external_gene_name)
mgi$H_external_gene_name <- gsub("[*]", "", mgi$H_external_gene_name)


# subset mouse models of human disease
# summary(factor(mgi$Gene_Catagory))
HMdg <- mgi[mgi$Gene_Catagory  == "MouseAndHuman",]
Mdg <- mgi[mgi$Gene_Catagory  == "MouseOnly",]
Hdg <- mgi[mgi$Gene_Catagory  == "HumanOnly",]

# subset catagories with n > 50
summary(factor(HMdg$Disease_cat))
df$HMDG_all <- 0
df$HMDG_all[df$H_external_gene_name %in% na.omit(unique(HMdg$H_external_gene_name))] <- 1
summary(factor(df$HMDG_all))

summary(factor(Hdg$Disease_cat))
df$HDG_all <- 0
df$HDG_all[df$H_external_gene_name %in% na.omit(unique(Hdg$H_external_gene_name))] <- 1
summary(factor(df$HDG_all))


### PLOT 
H_box.df <- df[,c("H_external_gene_name", 
                "H_fun_Z_percentile_0.001",
                     "HDG_all", "HMDG_all")]
H_box.df <- H_box.df[complete.cases(H_box.df),]
H_box.df <- H_box.df[!duplicated(H_box.df$H_external_gene_name),]
H_box.df <- H_box.df[H_box.df$HDG_all == 1,]
H.mm <- length(H_box.df$HMDG_all[H_box.df$HMDG_all == 1])
H.no.mm <- nrow(H_box.df) - H.mm
H_box.df$cat <- paste0("No mouse model\n(n=", H.no.mm, ")")
H_box.df$cat[H_box.df$HMDG_all == 1] <- paste0("Mouse model\n(n=", H.mm, ")")
summary(factor(H_box.df$cat))

my_comparisons <- list( c("Mouse model\n(n=296)", "No mouse model\n(n=1589)"))


H.box <- ggplot(H_box.df, aes(x=cat, y=H_fun_Z_percentile_0.001, fill=cat)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = my_comparisons,
    method = "wilcox.test", 
    label = "p.signif",
    paired = F,
    label.y = c(1.1)
    # size = 3,
  ) +
  ylab("Human funZ percentile") +
  xlab("") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1,0.5,1,1),"cm")) 
H.box

M_box.df <- df[,c("M_external_gene_name", 
                  "M_fun_Z_percentile_0.0001",
                  "HDG_all", "HMDG_all")]
M_box.df <- M_box.df[complete.cases(M_box.df),]
M_box.df <- M_box.df[!duplicated(M_box.df$M_external_gene_name),]
M_box.df <- M_box.df[M_box.df$HDG_all == 1,]
M.mm <- length(M_box.df$HMDG_all[M_box.df$HMDG_all == 1])
M.no.mm <- nrow(M_box.df) - M.mm
M_box.df$cat <- paste0("No mouse model\n(n=", M.no.mm, ")")
M_box.df$cat[M_box.df$HMDG_all == 1] <- paste0("Mouse model\n(n=", M.mm, ")")
summary(factor(M_box.df$cat))

my_comparisons <- list( c("Mouse model\n(n=300)", "No mouse model\n(n=1642)"))


M.box <- ggplot(M_box.df, aes(x=cat, y=M_fun_Z_percentile_0.0001, fill=cat)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = my_comparisons,
    method = "wilcox.test", 
    label = "p.signif",
    paired = F,
    label.y = c(1.1)
    # size = 3,
  ) +
  ylab("Mouse funZ percentile") +
  xlab("") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1,1,1,0.5),"cm")) 
M.box

Fig <- plot_grid(H.box, M.box, ncol = 2, nrow = 1)
Fig
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_4.jpg", Fig, ncol = 1, nrow = 1, base_height = 6, base_width = 9)




#########
# 
# ### PLOT 
# # one2one <- df[df$orthology_type == "ortholog_one2one",]
# box.df <- df[,c("H_external_gene_name", "M_fun_Z_percentile_0.0001",
#                   "H_fun_Z_percentile_0.001",
#                   "HDG_all", "HMDG_all")]
# box.df <- box.df[complete.cases(box.df),]
# box.df <- box.df[!duplicated(box.df$H_external_gene_name),]
# n.all <- length(unique(box.df$H_external_gene_name))
# n.HDG <- length(box.df$HDG_all[box.df$HDG_all == 1])
# n.HMDG <- length(box.df$HMDG_all[box.df$HMDG_all == 1])
# 
# box.df$variable <- "all"
# s1 <- subset(box.df, box.df$HDG_all == 1)
# s1$variable <- "HDG_all"
# s2 <- subset(box.df, box.df$HMDG_all == 1)
# s2$variable <- "HMDG_all"
# box.df <- rbind(box.df, s1, s2)
# box.df <- box.df[,c("H_external_gene_name", "H_fun_Z_percentile_0.001", "M_fun_Z_percentile_0.0001",
#                     "variable")]
# box.df <- melt(box.df, id.vars = c("H_external_gene_name", "variable"))
# box.df$variable.1 <- as.character(box.df$variable.1)
# box.df$variable.1[box.df$variable.1 == "H_fun_Z_percentile_0.001"] <- "Human genes"
# box.df$variable.1[box.df$variable.1 == "M_fun_Z_percentile_0.0001"] <- "Mouse orthologs"
# box.df$variable <- as.character(box.df$variable)
# box.df$variable[box.df$variable == "all"] <- "All genes\n(n = 14982)"
# box.df$variable[box.df$variable == "HDG_all"] <- "Human disease\ngenes\n(n = 1940)"
# box.df$variable[box.df$variable == "HMDG_all"] <- "Human disease\ngenes with\nmouse models\n(n = 858)"
# 
# 
# H_box.df <- subset(box.df, box.df$variable.1 == "Human genes")
# M_box.df <- subset(box.df, box.df$variable.1 == "Mouse orthologs")
# my_comparisons <- list( c("All genes\n(n = 14982)", "Human disease\ngenes\n(n = 1940)"),
#                         c("Human disease\ngenes\n(n = 1940)", "Human disease\ngenes with\nmouse models\n(n = 858)"))
# 
# H.box <- ggplot(H_box.df, aes(x=variable, y=value, fill=variable)) +
#   geom_boxplot() +
#   stat_compare_means(
#     aes(group = variable),
#     comparisons = my_comparisons,
#     method = "wilcox.test", 
#     label = "p.signif",
#     paired = F,
#     # size = 3,
#     label.y = c(1.1, 1.2, 1.3)
#   ) +
#   ylab("Human funZ percentile") +
#   xlab("") +
#   theme_classic() +
#   theme(legend.title=element_blank(),
#         legend.key.size = unit(1.5, "cm"),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.text=element_text(size=12),
#         legend.position="none",
#         plot.margin=unit(c(1,0.5,1,1),"cm")) 
# H.box
# 
# M.box <- ggplot(M_box.df, aes(x=variable, y=value, fill=variable)) +
#   geom_boxplot() +
#   stat_compare_means(
#     aes(group = variable),
#     comparisons = my_comparisons,
#     method = "wilcox.test", 
#     label = "p.signif",
#     paired = F,
#     # size = 3,
#     label.y = c(1.1, 1.2, 1.3)
#   ) +
#   ylab("Mouse funZ percentile") +
#   xlab("") +
#   theme_classic() +
#   theme(legend.title=element_blank(),
#         legend.key.size = unit(1.5, "cm"),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.text=element_text(size=12),
#         legend.position="none",
#         plot.margin=unit(c(1,1,1,0.5),"cm")) 
# M.box
# 
# Fig <- plot_grid(H.box, M.box, ncol = 2, nrow = 1)
# Fig
# save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_4.jpg", Fig, ncol = 1, nrow = 1, base_height = 8, base_width = 10)
# 

# out.box <- ggplot(box.df, aes(x=variable.1, y=value, fill=variable)) +
#   geom_boxplot() +
#   stat_compare_means(
#     aes(group = variable),
#     # comparisons = my_comparisons,
#     method = "wilcox.test", 
#     # label = "p.signif", 
#     paired = F,
#     hide.ns = T,
#     label.y = c(1.1, 1.2, 1.3)
#     # size = 8
#     ) +
#   ylab("funZ percentile") +
#   xlab("") +
#   theme_classic() +
#   theme(legend.title=element_blank(),
#         legend.key.size = unit(1.5, "cm"),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         legend.text=element_text(size=14),
#         # legend.position="top",
#         plot.margin=unit(c(1,1,1,1),"cm")) 
# out.box

# ggboxplot(box.df, x = "variable", y = "value",
#           color = "variable", 
#           facet.by = "variable.1", short.panel.labs = FALSE)

# wilcox.test(one2one$M_fun_Z_percentile_0.0001[one2one$HDG_all == 1], one2one$M_fun_Z_percentile_0.0001[one2one$HMDG_all == 1])
# x <- one2one[one2one$HDG_all == 1,]
# wilcox.test(x$M_fun_Z_percentile_0.0001[x$HMDG_all == 0], x$M_fun_Z_percentile_0.0001[x$HMDG_all == 1])

