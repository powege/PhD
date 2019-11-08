rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

### IMPORT

dt1 <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_annotation_alignment_by_sequence.csv")
dt2 <- fread("~/Dropbox/PhD/Data/Ensembl/Alignment/Formatted/HM_annotation_alignment_totals.csv")
summary <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/HM_annotation_coverage_summary.csv")

### FORMAT

# subset human annotations
summary <- summary[species == "Human"]
summary <- summary[,c("category", "genomic_percentage")]
summary <- rbind(summary, 
                 data.table(category = "Miscellaneous", 
                            genomic_percentage = (summary$genomic_percentage[summary$category == "TF binding"] + 
                                                  summary$genomic_percentage[summary$category == "Open chromatin"])))
summary <- summary[category != "TF binding" & category != "Open chromatin"]
summary$x_lab <- paste0(summary$category, "\n(", summary$genomic_percentage, "%)")

# rename annotations
dt1$category[dt1$category == "TF binding" | dt1$category == "Open chromatin"] <- "Miscellaneous"

# set factor order
order <- c("Exon - CDS",
           "Exon - UTR",
           "Exon - other",
           "Promoter",
           "Intron - proximal",
           "Enhancer - proximal",
           "Enhancer - distal",
           "CTCF binding",
           "Miscellaneous",
           "Intron - distal",
           "Unannotated")
order <- rev(order)
dt1$category <- factor(dt1$category, levels = as.character(order))
dt2$H_annotation <- factor(dt2$H_annotation, levels = as.character(order))

summary$category <- factor(summary$category, levels = as.character(order))
summary <- summary[order(category)]
lab_order <- summary$x_lab
summary$x_lab <- factor(summary$x_lab, levels = as.character(lab_order))

# format annotation conservation
dt_ann_con <- dt2[dt2$H_annotation == dt2$M_annotation,]
dt_ann_con$conserved_same_ann <- dt_ann_con$HM_ann_align_total / dt_ann_con$H_ann_total
dt_ann_con <- dt_ann_con[,c("H_annotation", "conserved_same_ann", "H_ann_align_frac")]
colnames(dt_ann_con) <- c("category", "total_aligned_same_ann", "total_aligned")
dt_ann_con$total_aligned_diff_ann <- dt_ann_con$total_aligned - dt_ann_con$total_aligned_same_ann
dt_ann_con <- melt(dt_ann_con, measure.vars = c("total_aligned_same_ann", "total_aligned_diff_ann"))
tmp <- aggregate(dt1[, "align_frac"], list(dt1$category), median)
colnames(tmp) <- c("category", "Q2")
dt_ann_con <- dt_ann_con[tmp, on = "category"]
dt_ann_con$variable <- factor(dt_ann_con$variable, levels = c("total_aligned_same_ann", "total_aligned_diff_ann"))
dt_ann_con$col <- ifelse(dt_ann_con$variable == "total_aligned_diff_ann", "black", as.character(dt_ann_con$category))
dt_ann_con$col <- factor(dt_ann_con$col, levels = c("black", as.character(order)))

dt_ann_con <- summary[dt_ann_con, on = "category"]

color_scale <- c("black",
                 brewer.pal(length(unique(dt_ann_con$col)) - 1, name = "Set3")
)
color_map <- setNames(color_scale, levels(unique(dt_ann_con$col)))


### PLOT A

pA <- ggplot(dt1, aes(x=category, y=align_frac, fill=category)) +
  geom_violin() +
  # geom_boxplot(width=0.1, outlier.shape = NA, coef = 0) +
  geom_point(data = dt_ann_con, aes(x=category, y=Q2), colour="black", shape=4, size=3, stroke=2) +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  ggtitle("A") +
  scale_fill_brewer(palette="Set3") +
  # geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pA

### PLOT B

pB <- ggplot(dt_ann_con, aes(x=category, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  ggtitle("B") +
  scale_fill_manual(values = color_map) +
  # geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pB


### OUTPUT

poutAB <- grid.arrange(pA, pB, nrow = 1, widths = c(6, 6))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_AB.jpg", plot = poutAB, height = 5, width = 12)

#####

### PLOT C

pC <- ggplot(dt_ann_con, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  # ggtitle("B") +
  scale_fill_manual(values = color_map) +
  # geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pC
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_C.jpg", plot = pC, height = 5, width = 6)
