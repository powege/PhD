rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(plyr)

### FUNCTION

alakazam <- function(dt, grouping){

colnames(dt) <- c("H_CHR", "H_POS", "H_REF", "H_KMER", "H_ANN", 
                  "M_CHR", "M_POS", "M_REF", "M_KMER", "M_ANN")

ann_short <- LETTERS[1:11]

n_SNV <- rep(NA, length(ann_short))
n_align <- rep(NA, length(ann_short))
n_align_ann <- rep(NA, length(ann_short))

for (i in 1:length(ann_short)){
  ann_sub <- dt[H_ANN %like% ann_short[i]]
  n_SNV[i] <- nrow(ann_sub)
  n_align[i] <- nrow(ann_sub[!is.na(M_POS)])
  n_align_ann[i] <- nrow(ann_sub[M_ANN %like% ann_short[i]])
}
 
total_align <- nrow(dt[!is.na(M_POS)]) / nrow(dt)

dt_out <- data.table(group = rep(grouping, length(ann_short)),
                     H_ANN = ann_short,
                     n_SNV = n_SNV,
                     n_align = n_align,
                     n_align_ann = n_align_ann,
                     frac_align = n_align/n_SNV,
                     frac_align_ann = n_align_ann/n_SNV,
                     total_align = rep(total_align, length(ann_short)))

dt_out$H_ANN <- mapvalues(dt_out$H_ANN, from=c("A","B","C","D","E","F","G","H","I","J", "K"), 
                       to=c("Exon - CDS",
                            "Exon - UTR",
                            "Exon - other",
                            "Promoter",
                            "Enhancer - proximal",
                            "Enhancer - distal",
                            "CTCF binding",
                            "Miscellaneous",
                            "Intron - proximal",
                            "Intron - distal",
                            "Unannotated"))

return(dt_out)
}

kadabra <- function(dt){
  
  dt$x_lab <- paste0(dt$H_ANN, "\n(n=", dt$n_SNV, ")")
  
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
  dt$H_ANN <- factor(dt$H_ANN, levels = as.character(order))
  dt <- dt[order(H_ANN)]
  lab_order <- dt$x_lab
  dt$x_lab <- factor(dt$x_lab, levels = as.character(lab_order))
  
  dt$frac_align_diff_ann <- dt$frac_align - dt$frac_align_ann
  
  dt_melt <- melt(dt, measure.vars = c("frac_align_diff_ann", "frac_align_ann"))
  dt_melt$variable <- factor(dt_melt$variable, levels = c("frac_align_ann", "frac_align_diff_ann"))
  dt_melt$col <- ifelse(dt_melt$variable == "frac_align_diff_ann", "black", as.character(dt_melt$H_ANN))
  dt_melt$col <- factor(dt_melt$col, levels = c("black", as.character(order)))
  
  return(dt_melt)
}

### IMPORT

nonDDC <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_nonDDC_mouse_mapping.csv")
DDC <- fread("~/Dropbox/PhD/Data/GWAS/formatted/GWAS_DDC_mouse_mapping.csv")
CV <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_mouse_mapping.csv")

### FORMAT

nonDDC_out <- alakazam(dt = nonDDC, grouping = "GWAS nonDDC")
DDC_out <- alakazam(dt = DDC, grouping = "GWAS DDC")
CV_out <- alakazam(dt = CV, grouping = "ClinVar")

nonDDC_plot <- kadabra(dt = nonDDC_out)
DDC_plot <- kadabra(dt = DDC_out)
CV_plot <- kadabra(dt = CV_out)


### PLOT 

color_scale <- c("black",
                 brewer.pal(length(unique(DDC_plot$col)) - 1, name = "Set3")
)
color_map <- setNames(color_scale, levels(unique(DDC_plot$col)))

p_DDC <- ggplot(DDC_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  ggtitle("GWAS") +
  scale_fill_manual(values = color_map) +
  # geom_hline(yintercept=DDC_plot$total_align[1], linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_DDC
ggsave("~/Dropbox/PhD/Data/Figures/Figure_DDC.jpg", plot = p_DDC, height = 5, width = 6)


p_nonDDC <- ggplot(nonDDC_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  # ggtitle("B") +
  scale_fill_manual(values = color_map) +
  # geom_hline(yintercept=nonDDC_plot$total_align[1], linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_nonDDC
ggsave("~/Dropbox/PhD/Data/Figures/Figure_nonDDC.jpg", plot = p_nonDDC, height = 5, width = 6)

p_CV <- ggplot(CV_plot, aes(x=x_lab, y=value, fill=col)) +
  geom_bar(stat="identity", colour = "black") +
  coord_flip() +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  ggtitle("ClinVar") +
  scale_fill_manual(values = color_map) +
  # geom_hline(yintercept=CV_plot$total_align[1], linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p_CV
ggsave("~/Dropbox/PhD/Data/Figures/Figure_ClinVar.jpg", plot = p_CV, height = 5, width = 6)


pout <- grid.arrange(p_DDC, p_CV, nrow = 1, widths = c(6, 6))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_ClinVar_GWAS.jpg", plot = pout, height = 5, width = 12)












 