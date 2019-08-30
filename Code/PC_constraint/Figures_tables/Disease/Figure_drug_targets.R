rm(list = ls())
graphics.off()

library(data.table)
library(plotROC)
library(pROC)
library("cowplot")

##############
### IMPORT RAW
##############

drug <- fread("~/Dropbox/PhD/Data/Drug_Targets/Drug_targets_list_for_George_P.csv", fill = T, header = F)
funz <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
rvis <- fread("~/Dropbox/PhD/Data/Constraint/RVIS_Petrovski2013_ESP.csv")
gnomad <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt")

##############################
### FORMAT constraint metrics
#############################

gnomad <- gnomad[,c("gene", "mis_z", "pLI", "oe_lof_upper")]
colnames(gnomad) <- c("H_external_gene_name", "mis_z", "pLI", "oe_lof_upper")
rvis <- rvis[,c("HGNC gene", "Residual Variation Intolerance Score")]
colnames(rvis) <- c("H_external_gene_name", "RVIS")
funz <- funz[,c("H_external_gene_name", "M_external_gene_name", "M_fun_Z", "H_fun_Z", "orthology_type")]

# merge
metrics <- gnomad[rvis, on = "H_external_gene_name"]
metrics <- metrics[funz, on = "H_external_gene_name"]
metrics <- metrics[!duplicated(metrics),]
# metrics <- metrics[complete.cases(metrics),]

# calculate percentiles
percentile <- ecdf(metrics$mis_z[!duplicated(metrics$H_external_gene_name)])
metrics$mis_z_percentile <- percentile(metrics$mis_z)
percentile <- ecdf(metrics$pLI[!duplicated(metrics$H_external_gene_name)])
metrics$pLI_percentile <- percentile(metrics$pLI)
percentile <- ecdf(metrics$oe_lof_upper[!duplicated(metrics$H_external_gene_name)])
metrics$oe_lof_upper_percentile <- percentile(metrics$oe_lof_upper)
percentile <- ecdf(metrics$RVIS[!duplicated(metrics$H_external_gene_name)])
metrics$RVIS_percentile <- percentile(metrics$RVIS)
percentile <- ecdf(metrics$M_fun_Z[!duplicated(metrics$M_external_gene_name)])
metrics$M_fun_Z_percentile <- percentile(metrics$M_fun_Z)
percentile <- ecdf(metrics$H_fun_Z[!duplicated(metrics$H_external_gene_name)])
metrics$H_fun_Z_percentile <- percentile(metrics$H_fun_Z)

# one_to_one orthologues
# metrics <- subset(metrics, metrics$orthology_type == "ortholog_one2one")

######################
### FORMAT gene lists
#####################

## Identify genes in each gene list
drug_t <- unlist(strsplit(drug$V1, ","))
length(drug_t)

## Identify genes not in drug list
drug_t_neg <- unique(metrics$H_external_gene_name[which(!metrics$H_external_gene_name %in% drug_t)])
length(drug_t_neg)

gene_list_dt <- data.table(H_external_gene_name = c(drug_t, drug_t_neg),
                           sig = c(rep(1, length(drug_t)), rep(0, length(drug_t_neg))))


##############
### WORKING DT
##############

dt <- gene_list_dt[metrics, on = "H_external_gene_name"]
dt_long <- dt[,c("H_external_gene_name", "sig", "mis_z_percentile",      
                 "pLI_percentile", "oe_lof_upper_percentile", "RVIS_percentile",      
                 "M_fun_Z_percentile", "H_fun_Z_percentile")] 
dt_long <- melt(dt_long, 
                id.vars = c("H_external_gene_name","sig"), 
                measure.vars = c("mis_z_percentile", "pLI_percentile", 
                                 "oe_lof_upper_percentile", "RVIS_percentile", 
                                 "M_fun_Z_percentile", "H_fun_Z_percentile"))
dt_long <- dt_long[complete.cases(dt_long),]
dt_long$sig <- as.factor(dt_long$sig)

###############
### STAT TEST 
############

wilcox.test(dt_long$value[dt_long$variable == "M_fun_Z_percentile" & dt_long$sig == 1],
            dt_long$value[dt_long$variable == "M_fun_Z_percentile" & dt_long$sig == 0])
wilcox.test(dt_long$value[dt_long$variable == "H_fun_Z_percentile" & dt_long$sig == 1],
            dt_long$value[dt_long$variable == "H_fun_Z_percentile" & dt_long$sig == 0])
wilcox.test(dt_long$value[dt_long$variable == "RVIS_percentile" & dt_long$sig == 1],
            dt_long$value[dt_long$variable == "RVIS_percentile" & dt_long$sig == 0])
wilcox.test(dt_long$value[dt_long$variable == "oe_lof_upper_percentile" & dt_long$sig == 1],
            dt_long$value[dt_long$variable == "oe_lof_upper_percentile" & dt_long$sig == 0])
wilcox.test(dt_long$value[dt_long$variable == "pLI_percentile" & dt_long$sig == 1],
            dt_long$value[dt_long$variable == "pLI_percentile" & dt_long$sig == 0])
wilcox.test(dt_long$value[dt_long$variable == "mis_z_percentile" & dt_long$sig == 1],
            dt_long$value[dt_long$variable == "mis_z_percentile" & dt_long$sig == 0])

##################
### FIGURE boxplot
##################

fig1 <- ggplot(data = dt_long, aes(x=variable, y=value, fill = sig)) +
  # geom_boxplot() +
  # geom_rect(
  #   aes(xmin = dt$xmin,
  #       xmax = dt$xmax,
  #       ymin = -Inf,
  #       ymax = Inf,
  #       fill = dt$rect_type),
  #   alpha = 0.01) +
  # stat_compare_means(label = "p.signif",
  #                    method = "t.test",
  #                    symnum.args = symnum.args,
#                    size = 6,
#                    hide.ns = T) +
geom_boxplot() +
  # stat_compare_means(
  #   aes(group = sig),
  #   # comparisons = comparisons,
  #   method = "wilcox.test", 
  #   label = "p.signif"
  #   # label.y = c(1.15, 1.25, 1.05),
  #   # paired = F
  #   # size = 3,
  # ) +
  # scale_fill_manual(values=c("skyblue")) +
  # scale_fill_brewer(palette="Set2") +
  coord_flip() +
  xlab("") +
  ylab("Constraint percentile rank") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = "none",
        text = element_text(size = 14)) 
fig1

