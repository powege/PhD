rm(list=ls())
graphics.off()

library(corrplot)
library(tidyr)
library(RColorBrewer)
library("cowplot")
library(data.table)

### IMPORT

funZ <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
M_funZ_AS <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_MGP_all_strains_v2.csv")
H_funZ_0001 <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF0001_v2.csv")
H_funZ_0005 <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_1KGP_MAF0005_v2.csv")
H_funZ_gAD <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/H_funZ_gnomAD_MAF001_v2.csv")


### FORMAT

funZ <- funZ[,c("H_external_gene_name", "M_external_gene_name", 
                "M_fun_Z", "H_fun_Z", 
                "orthology_type")]
M_funZ_AS <- M_funZ_AS[,c("external_gene_name", "fun_Z")]
colnames(M_funZ_AS) <- c("M_external_gene_name", "M_fun_Z_all_strains")

H_funZ_0001 <- H_funZ_0001[,c("external_gene_name", "fun_Z")]
colnames(H_funZ_0001) <- c("H_external_gene_name", "H_fun_Z_MAF0001")

H_funZ_0005 <- H_funZ_0005[,c("external_gene_name", "fun_Z")]
colnames(H_funZ_0005) <- c("H_external_gene_name", "H_fun_Z_MAF0005")

H_funZ_gAD <- H_funZ_gAD[,c("external_gene_name", "fun_Z")]
colnames(H_funZ_gAD) <- c("H_external_gene_name", "H_fun_Z_gnomAD")

# subset one to one orthologues
funZ <- subset(funZ, funZ$orthology_type == "ortholog_one2one")

# merge
dt_all2 <- H_funZ_gAD[H_funZ_0005, on = "H_external_gene_name"]
dt_all2 <- dt_all2[H_funZ_0001, on = "H_external_gene_name"]
dt_all2 <- dt_all2[funZ, on = "H_external_gene_name"]
dt_all2 <- dt_all2[M_funZ_AS, on = "M_external_gene_name"]
dt_all2 <- dt_all2[!duplicated(dt_all2),]
dt_all2 <- dt_all2[complete.cases(dt_all2),]

# correlation matrix
dt_plot2 <- cor(dt_all2[,c("M_fun_Z",
                           "M_fun_Z_all_strains",
                           "H_fun_Z",
                           "H_fun_Z_MAF0005",
                           "H_fun_Z_MAF0001",
                           "H_fun_Z_gnomAD")], method = "spearman")
dt_plot2 <- as.matrix(dt_plot2)
colnames(dt_plot2) <- c("Mouse funZ no SPRET/EiJ",
                        "Mouse funZ all strains",
                        "Human funZ 1KGP MAF>001",
                        "Human funZ 1KGP MAF>0005",
                        "Human funZ 1KGP MAF>0001",
                        "Human funZ gnomAD MAF>001")
rownames(dt_plot2) <- c("Mouse funZ no SPRET/EiJ",
                        "Mouse funZ all strains",
                        "Human funZ 1KGP MAF>001",
                        "Human funZ 1KGP MAF>0005",
                        "Human funZ 1KGP MAF>0001",
                        "Human funZ gnomAD MAF>001")

### PLOT

pdf(file = "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_constraint_cor_matrix_supplement.pdf", width = 6, height = 6)
corrplot(dt_plot2,
         method = "color",
         # type = "upper",
         # cl.lim = c(0, 1),
         col = brewer.pal(n = 10, name = "RdYlBu"),
         mar = c(0,0,0,3),
         number.cex = 0.9,
         number.digits = 2,
         addCoef.col = "black",
         tl.col = "black")
dev.off()


########

# # jpeg(file = "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_constraint_cor.jpg", 
# #      width = 12, height = 6, units = "in", res = 1000)
# pdf(file = "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_constraint_cor_matrix.pdf", width = 6, height = 6)
# corrplot(dt_plot,
#          method = "color",
#          # type = "upper",
#          # cl.lim = c(0, 1),
#          col = brewer.pal(n = 10, name = "RdYlBu"),
#          mar = c(0,0,0,3),
#          number.cex = 0.9,
#          number.digits = 2,
#          addCoef.col = "black",
#          tl.col = "black")
# dev.off()
