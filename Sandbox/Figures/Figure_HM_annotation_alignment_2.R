rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

### IMPORT

dt1 <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_alignment_by_human_annotation.csv")
dt2 <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_alignment_annotation_summary.csv")
dt3H <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/H_annotation_overlap_summary.csv")
dt3M <- fread("~/Dropbox/PhD/Data/Ensembl/Annotation/M_annotation_overlap_summary.csv")

### FORMAT

# set factor order
order <- c("Exon - CDS",
           "Exon - UTR",
           "Exon - other",
           "Promoter",
           "Intron - proximal",
           "Enhancer - proximal",
           "Enhancer - distil",
           "CTCF binding",
           "Open chromatin",
           "TF binding",
           "Intron - distil",
           "Unannotated")
order <- rev(order)
dt1$category <- factor(dt1$category, levels = as.character(order))
dt2$H_annotation <- factor(dt2$H_annotation, levels = as.character(order))
dt3H$annotation_primary <- factor(dt3H$annotation_primary, levels = as.character(order))
dt3H$annotation_secondary <- factor(dt3H$annotation_secondary, levels = as.character(order))
dt3M$annotation_primary <- factor(dt3M$annotation_primary, levels = as.character(order))
dt3M$annotation_secondary <- factor(dt3M$annotation_secondary, levels = as.character(order))

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

color_scale <- c("black",
                 brewer.pal(length(unique(dt_ann_con$col)) - 1, name = "Set3")
)
color_map <- setNames(color_scale, levels(dt_ann_con$col))


### FORMAT dt3 

dt3H_wide <- split(dt3H , f = dt3H$chromosome)
dt3H_wide_p_total <- lapply(dt3H_wide, "[[", 4)
dt3H_wide_p_total <- do.call("cbind", dt3H_wide_p_total)
dt3H_wide_s_total <- lapply(dt3H_wide, "[[", 5)
dt3H_wide_s_total <- do.call("cbind", dt3H_wide_s_total)

dt_pC <- data.table(annotation_primary = dt3H_wide[[1]]$annotation_primary,
                    annotation_secondary = dt3H_wide[[1]]$annotation_secondary,
                    primary_total_bp = apply(dt3H_wide_p_total, 1, sum),
                    secondary_total_bp = apply(dt3H_wide_s_total, 1, sum))
dt_pC$fraction <- dt_pC$secondary_total_bp/dt_pC$primary_total_bp
dt_pC <- subset(dt_pC, dt_pC$annotation_primary != dt_pC$annotation_secondary)

dt3M_wide <- split(dt3M , f = dt3M$chromosome)
dt3M_wide_p_total <- lapply(dt3M_wide, "[[", 4)
dt3M_wide_p_total <- do.call("cbind", dt3M_wide_p_total)
dt3M_wide_s_total <- lapply(dt3M_wide, "[[", 5)
dt3M_wide_s_total <- do.call("cbind", dt3M_wide_s_total)

dt_pD <- data.table(annotation_primary = dt3M_wide[[1]]$annotation_primary,
                    annotation_secondary = dt3M_wide[[1]]$annotation_secondary,
                    primary_total_bp = apply(dt3M_wide_p_total, 1, sum),
                    secondary_total_bp = apply(dt3M_wide_s_total, 1, sum))
dt_pD$fraction <- dt_pD$secondary_total_bp/dt_pD$primary_total_bp
dt_pD <- subset(dt_pD, dt_pD$annotation_primary != dt_pD$annotation_secondary)

### Whole genome alignment
GW_Halignment <- 0.2926328

### Percent of genome
# set human chr lengths (GRCh38) from 
h_chr_len <- c(248956422,	242193529,198295559,190214555,181538259,170805979,159345973,145138636,
               138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
               83257441,80373285,58617616,64444167,46709983,50818468)
tmp_sub <- dt_pC[,c("annotation_primary", "primary_total_bp")]
tmp_sub <- tmp_sub[!duplicated(tmp_sub),]
tmp_sub$coverage <- round( ((tmp_sub$primary_total_bp / sum(h_chr_len)) * 100), digits = 1)

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
  geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
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
  geom_hline(yintercept=GW_Halignment, linetype="dashed", color = "blue", size=1) +
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

### PLOT C 

pC <- ggplot() +
  geom_bar(data=dt_pC, aes(x=annotation_primary, y=fraction, fill=annotation_secondary), colour = "black", stat="identity") +
  xlab("Genomic annotation") +
  ylab("Intraspecific overlap (fraction)") +
  ggtitle("Human") +
  # labs(fill = "") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     # trans = "reverse",
                     limits = c(0, 0.85)) +
  coord_flip() +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    # legend.title = element_text(hjust = 0.5),
    # legend.key.size = unit(2, 'lines'),
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.05),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pC

pD <- ggplot() +
  geom_bar(data=dt_pD, aes(x=annotation_primary, y=fraction, fill=annotation_secondary), colour = "black", stat="identity") +
  xlab("Genomic annotation") +
  ylab("Intraspecific overlap (fraction)") +
  ggtitle("Mouse") +
  # labs(fill = "") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     # trans = "reverse",
                     limits = c(0, 0.85)) +
  coord_flip() +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    # legend.title = element_text(hjust = 0.5),
    # legend.key.size = unit(2, 'lines'),
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.05),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pD


### OUTPUT

# pout <- grid.arrange(pA, pB, nrow = 2)
# ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_tmp.jpg", plot = pout, height = 12, width = 6)
poutAB <- grid.arrange(pA, pB, nrow = 1, widths = c(6, 6))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_AB.jpg", plot = poutAB, height = 5, width = 12)

poutCD <- grid.arrange(pC, pD, nrow = 1, widths = c(5, 7))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_CD.jpg", plot = poutCD, height = 5, width = 12)

pout <- grid.arrange(poutAB, poutCD, nrow = 2)
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_tmp.jpg", plot = pout, height = 10, width = 12)


