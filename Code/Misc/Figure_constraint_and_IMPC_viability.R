rm(list = ls())
graphics.off()

library(ggplot2)
library(data.table)
library(ggpubr)



### IMPORT 

load("Dropbox/allLinesWithPValues.RData")
funz <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
gnomad <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt")

### FORMAT

gnomad <- gnomad[,c("gene", "pLI", "oe_lof_upper")]
colnames(gnomad) <- c("H_external_gene_name", "pLI", "oe_lof_upper")

percentile <- ecdf(gnomad$pLI[!duplicated(gnomad$H_external_gene_name)])
gnomad$pLI_percentile <- percentile(gnomad$pLI)
percentile <- ecdf(gnomad$oe_lof_upper[!duplicated(gnomad$H_external_gene_name)])
gnomad$oe_lof_upper_percentile <- percentile(gnomad$oe_lof_upper)

lines <- lines[,c("gene_symbol", "outcomeFinal")]
colnames(lines) <- c("M_external_gene_name", "outcome")

funz <- funz[,c("H_external_gene_name", "M_external_gene_name", 
                "M_fun_Z", 
                "M_OE_ratio")]

percentile <- ecdf(funz$M_fun_Z[!duplicated(funz$M_external_gene_name)])
funz$fun_Z_percentile <- percentile(funz$M_fun_Z)
percentile <- ecdf(funz$M_OE_ratio[!duplicated(funz$M_external_gene_name)])
funz$OE_ratio_percentile <- percentile(funz$M_OE_ratio)

funz$OE_ratio_percentile <- 1-funz$OE_ratio_percentile
gnomad$oe_lof_upper_percentile <- 1-gnomad$oe_lof_upper_percentile


plot_dt <- funz[lines, on = "M_external_gene_name"]
plot_dt <- gnomad[plot_dt, on = "H_external_gene_name"]

plot_dt <- plot_dt[,c("H_external_gene_name", "M_external_gene_name", "outcome", 
                      # "OE_ratio_percentile",
                      "fun_Z_percentile",
                      # "oe_lof_upper_percentile",
                      "pLI_percentile")]
plot_dt <- plot_dt[complete.cases(plot_dt),]
plot_dt <- melt(plot_dt, id.vars = c("H_external_gene_name", "M_external_gene_name", "outcome"))

plot_dt$variable <- as.character(plot_dt$variable)
plot_dt$variable[plot_dt$variable == "fun_Z_percentile"] <- "Mouse funZ"
plot_dt$variable[plot_dt$variable == "pLI_percentile"] <- "Human pLI"
plot_dt$outcome <- as.character(plot_dt$outcome)
plot_dt$outcome[plot_dt$outcome == "Severe Subviable"] <- "Severe\nSubviable"
plot_dt$outcome[plot_dt$outcome == "Statistical Subviable"] <- "Statistical\nSubviable"

plot_dt_pLI <- subset(plot_dt, plot_dt$variable == "Human pLI")
plot_dt_funZ <- subset(plot_dt, plot_dt$variable == "Mouse funZ")


### STAT TEST 

wilcox.test(plot_dt$value[plot_dt$outcome == "Lethal" & plot_dt$variable == "Mouse funZ"],
            plot_dt$value[plot_dt$outcome == "Severe\nSubviable" & plot_dt$variable == "Mouse funZ"])
wilcox.test(plot_dt$value[plot_dt$outcome == "Lethal" & plot_dt$variable == "Human pLI"],
            plot_dt$value[plot_dt$outcome == "Severe\nSubviable" & plot_dt$variable == "Human pLI"])

wilcox.test(plot_dt$value[plot_dt$outcome == "Severe\nSubviable" & plot_dt$variable == "Mouse funZ"],
            plot_dt$value[plot_dt$outcome == "Statistical\nSubviable" & plot_dt$variable == "Mouse funZ"])
wilcox.test(plot_dt$value[plot_dt$outcome == "Severe\nSubviable" & plot_dt$variable == "Human pLI"],
            plot_dt$value[plot_dt$outcome == "Statistical\nSubviable" & plot_dt$variable == "Human pLI"])

wilcox.test(plot_dt$value[plot_dt$outcome == "Statistical\nSubviable" & plot_dt$variable == "Mouse funZ"],
            plot_dt$value[plot_dt$outcome == "Viable" & plot_dt$variable == "Mouse funZ"])
wilcox.test(plot_dt$value[plot_dt$outcome == "Statistical\nSubviable" & plot_dt$variable == "Human pLI"],
            plot_dt$value[plot_dt$outcome == "Viable" & plot_dt$variable == "Human pLI"])


wilcox.test(plot_dt_funZ$value[plot_dt_funZ$outcome == "Statistical\nSubviable" | plot_dt_funZ$outcome == "Severe\nSubviable"],
            plot_dt_funZ$value[plot_dt_funZ$outcome == "Viable"])$p.val
wilcox.test(plot_dt_pLI$value[plot_dt_pLI$outcome == "Statistical\nSubviable" | plot_dt_pLI$outcome == "Severe\nSubviable"],
            plot_dt_pLI$value[plot_dt_pLI$outcome == "Viable"])$p.val

### PLOT

table(plot_dt_funZ$outcome)
M_comparisons <- list( c("Severe\nSubviable", "Statistical\nSubviable"))
                       # , 
                       # c("Statistical\nSubviable", "Viable"))
M_box <- ggplot(plot_dt_funZ, aes(x=outcome, y=value, fill = outcome)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = M_comparisons,
    method = "wilcox.test", 
    label = "p.format",
    label.y = c(1.15, 1.25, 1.05),
    paired = F
    # size = 3,
  ) +
  ylab("Mouse funZ percentile rank") +
  xlab("") +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_brewer(palette="Set2") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"cm")) 
M_box
ggsave(filename = "~/Dropbox/PhD/Data/Figures/Figure_IMPC_viability_funZ.jpg", plot = M_box, height = 6, width = 6)


table(plot_dt_pLI$outcome)
H_comparisons <- list( c("Severe\nSubviable", "Statistical\nSubviable"))
                       # ,
                       # c("Statistical\nSubviable", "Viable"))
H_box <- ggplot(plot_dt_pLI, aes(x=outcome, y=value, fill = outcome)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = H_comparisons,
    method = "wilcox.test", 
    label = "p.format",
    label.y = c(1.15, 1.25, 1.05),
    paired = F
    # size = 3,
  ) +
  ylab("Human pLI percentile rank") +
  xlab("") +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_brewer(palette="Set2") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"cm")) 
H_box
ggsave(filename = "~/Dropbox/PhD/Data/Figures/Figure_IMPC_viability_pLI.jpg", plot = H_box, height = 6, width = 6)

Fig <- plot_grid(M_box, H_box, ncol = 2, nrow = 1, labels = c("Mouse", "Human"), 
                 label_size = 20, hjust = -0.1, vjust = 2)
Fig
save_plot("~/Dropbox/PhD/Data/Figures/Figure_IMPC_viability_pLI_funZ.jpg", Fig, ncol = 1, nrow = 1, base_height = 5.5, base_width = 11)




box <- ggplot(plot_dt, aes(x=outcome, y=value, fill=variable)) +
  geom_boxplot() +
  ylab("Constraint percentile rank") +
  xlab("") +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_brewer(palette="Set2") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="top",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin=unit(c(1,1,1,1),"cm")) 
box
ggsave(filename = "~/Dropbox/PhD/Data/Figures/Figure_IMPC_viability.jpg", plot = box, height = 6, width = 6)




 #######

# ggplot(plot_dt, aes(x=outcome, y=fun_Z_IMPC_percentile, fill=outcome)) +
#     geom_boxplot()
# 
# ggplot(plot_dt, aes(x=outcome, y=fun_Z_percentile, fill=outcome)) +
#   geom_boxplot()
# 
# ggplot(plot_dt, aes(x=outcome, y=OE_ratio_IMPC_percentile, fill=outcome)) +
#   geom_boxplot()
# 
# ggplot(plot_dt, aes(x=outcome, y=OE_ratio_percentile, fill=outcome)) +
#   geom_boxplot()



