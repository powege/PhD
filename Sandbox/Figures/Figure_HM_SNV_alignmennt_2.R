rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)


### IMPORT
cv_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_z_scores_v2.csv")
gwas_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_z_scores_v2.csv")


### FORMAT

# add category
cv_z$SNV_category <- "ClinVar"
gwas_z$SNV_category <- "GWAS"

# subset plotting data tables
dt <- rbind(cv_z, gwas_z)
dt_p1 <- dt[, c("H_annotation", "SNV_category", "N_SNV_obs", "F_align_z")]
dt_p2 <- dt[, c("H_annotation", "SNV_category", "N_align_obs", "F_conserved_z")]
dt_p3 <- dt[, c("H_annotation", "SNV_category", "N_align_obs", "F_align_ann_z")]

# set factor order 
order_H <- c("Unannotated", "Intron", "TF binding", "Open chromatin", "Enhancer", "Promoter flanking",
             "Promoter", "Exon - non-coding", "Exon - UTR")
dt_p1$H_annotation  <- factor(dt_p1$H_annotation, levels = as.character(order_H))
dt_p2$H_annotation  <- factor(dt_p2$H_annotation, levels = as.character(order_H))
dt_p3$H_annotation  <- factor(dt_p3$H_annotation, levels = as.character(order_H))

# remove unannotated and intron
dt_p1 <- subset(dt_p1, dt_p1$H_annotation != "Intron" & dt_p1$H_annotation != "Unannotated")
dt_p2 <- subset(dt_p2, dt_p2$H_annotation != "Intron" & dt_p2$H_annotation != "Unannotated")
dt_p3 <- subset(dt_p3, dt_p3$H_annotation != "Intron" & dt_p3$H_annotation != "Unannotated")


### PLOT

p1 <- ggplot(dt_p1, aes(x=H_annotation, y=F_align_z, fill=SNV_category)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Human genomic annotation") +
  ylab("z-score") +
  ggtitle("Orthologous base in mouse") +
  coord_flip() +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20, 30),
                     limits = c(-20, 32)) +
  geom_hline(yintercept=3, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=-3, linetype="dashed", color = "black", size=1) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    # legend.position="bottom",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.65),
    # legend.box.background = element_rect(colour = "black"),
    # legend.spacing.x = unit(1.0, 'cm'),
    # legend.text = element_text(margin = margin(t = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) 
  # + guides(fill=guide_legend(
  #   keywidth=0.4,
  #   keyheight=0.4,
  #   default.unit="inch")
  # )
p1

p2 <- ggplot(dt_p2, aes(x=H_annotation, y=F_conserved_z, fill=SNV_category)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Human genomic annotation") +
  ylab("z-score") +
  ggtitle("Conserved base in mouse") +
  coord_flip() +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20, 30),
                     limits = c(-20, 32)) +
  geom_hline(yintercept=3, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=-3, linetype="dashed", color = "black", size=1) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    # legend.position="bottom",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.65),
    # legend.box.background = element_rect(colour = "black"),
    # legend.spacing.x = unit(1.0, 'cm'),
    # legend.text = element_text(margin = margin(t = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) 
p2

p3 <- ggplot(dt_p3, aes(x=H_annotation, y=F_align_ann_z, fill=SNV_category)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  xlab("Human genomic annotation") +
  ylab("z-score") +
  ggtitle("Conserved annotation in mouse") +
  coord_flip() +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20, 30),
                     limits = c(-20, 32)) +
  geom_hline(yintercept=3, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=-3, linetype="dashed", color = "black", size=1) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    # legend.position="bottom",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.65),
    # legend.box.background = element_rect(colour = "black"),
    # legend.spacing.x = unit(1.0, 'cm'),
    # legend.text = element_text(margin = margin(t = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) 
p3

pout <- grid.arrange(p1, p2, p3, nrow = 3)
ggsave("~/Dropbox/PhD/Data/Figures/Figure_SNV_Z_tmp2.jpg", plot = pout, height = 12, width = 6)


