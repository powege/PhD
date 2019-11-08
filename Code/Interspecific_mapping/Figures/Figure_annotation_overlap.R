rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

### IMPORT

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
dt3H$annotation_primary <- factor(dt3H$annotation_primary, levels = as.character(order))
dt3H$annotation_secondary <- factor(dt3H$annotation_secondary, levels = as.character(order))
dt3M$annotation_primary <- factor(dt3M$annotation_primary, levels = as.character(order))
dt3M$annotation_secondary <- factor(dt3M$annotation_secondary, levels = as.character(order))

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


### PLOT D

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


### PLOT A

dt_pA <- dt_pC
dt_pA$fraction <- round(dt_pA$fraction, digits = 2)
pA <- ggplot(dt_pA, aes(annotation_primary, annotation_secondary, fill = fraction))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red3", high = "blue3", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Overlap") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(vjust = 1, 
                                   size = 14, hjust = 1)) +
  coord_fixed() +
  ggtitle("Human") +
  xlab("Genomic annotation") +
  ylab("Intraspecific overlap (fraction)") +
  geom_text(aes(annotation_primary, annotation_secondary, label = fraction, fontface=2), color = "black", size = 5) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(1.3, 0.0),
    # legend.direction = "vertical"
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size = 14)
  )+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 12,
                               title.position = "top", title.hjust = 0.5))
pA


### PLOT B

dt_pB <- dt_pD
dt_pB$fraction <- round(dt_pB$fraction, digits = 2)
pB <- ggplot(dt_pB, aes(annotation_primary, annotation_secondary, fill = fraction))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red3", high = "blue3", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Overlap") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(vjust = 1, 
                                   size = 14, hjust = 1)) +
  coord_fixed() +
  ggtitle("Mouse") +
  xlab("Genomic annotation") +
  ylab("Intraspecific overlap (fraction)") +
  geom_text(aes(annotation_primary, annotation_secondary, label = fraction, fontface=2), color = "black", size = 5) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    # axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(1.3, 0.0),
    # legend.direction = "vertical"
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size = 14)
  )+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 12,
                               title.position = "top", title.hjust = 0.5))
pB


### OUTPUT

poutCD <- grid.arrange(pC, pD, nrow = 1, widths = c(5, 7))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_CD.jpg", plot = poutCD, height = 5, width = 12)

# poutAB <- grid.arrange(pA, pB, nrow = 2, ncol = 1)
# ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_overlap_AB.jpg", plot = poutAB, height = 20, width = 10)
poutAB <- grid.arrange(pA, pB, nrow = 1, ncol = 2)
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp_overlap_AB.jpg", plot = poutAB, height = 10, width = 20)


####

