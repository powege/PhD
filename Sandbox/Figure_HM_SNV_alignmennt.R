rm(list = ls())
graphics.off()

### ADD p-values

library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)

### FUNCTIONS

format_fun <- function(gwas_ddc_z){

# subset cols
p_ddc <- gwas_ddc_z[,c("H_annotation", "N_SNV_obs", "F_align_z", "F_conserved_z")]
# set x labs
p_ddc$x_lab <- paste0(p_ddc$H_annotation, "\n(n=", p_ddc$N_SNV_obs, ")")
# set factor order 
order_H <- unique(p_ddc$x_lab)
p_ddc$x_lab <- factor(p_ddc$x_lab, levels = as.character(rev(order_H)))

# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 11 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    x_lab = levels(p_ddc$x_lab),
                    rect_type = c("a", "c", "a", "c", "a", "c", "a", "c", "a", "c", "a"))
rects <- rbind(rects)
p_ddc <- merge(p_ddc, rects)
p_ddc <- melt(p_ddc, id.vars = c("x_lab", "H_annotation",  "N_SNV_obs", "xmin", "xmax", "rect_type"))

return(p_ddc)
}


plot_fun <- function(p_ddc, title){

pA <- ggplot(data = p_ddc, aes(x=x_lab, y=value, colour = variable)) + 
  geom_point(position=position_dodge(width=0.8), shape=4, size=3, stroke=2) +
  geom_bar(aes(x=x_lab, y=value, fill = variable), stat = "identity", position=position_dodge(width=0.8), width = 0.1) +
  geom_rect(
    aes(xmin = p_ddc$xmin,
        xmax = p_ddc$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = p_ddc$rect_type),
    color = NA,
    alpha = 0.5,
    show.legend = F) +
  geom_point(position=position_dodge(width=0.8), shape=4, size=3, stroke=2) +
  geom_bar(aes(x=x_lab, y=value, fill = variable), stat = "identity", position=position_dodge(width=0.8), width = 0.1) +
  scale_fill_manual(values = c("grey", "white", "blue", "red")) +
  scale_colour_manual(values = c("blue", "red")) +
  xlab("Human genomic annotation") +
  ylab("Z-score") +
  ggtitle(title) +
  geom_hline(yintercept=3, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=0, color = "black", size=0.8) +
  geom_hline(yintercept=-3, linetype="dashed", color = "black", size=1) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.title = element_blank(),
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.05),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
pA
return(pA)
}

### IMPORT
gwas_ddc_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_DDC_z_scores.csv")
gwas_non_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_nonDDC_z_scores.csv")

### FORMAT
ddc_dt <- format_fun(gwas_ddc_z)
non_dt <- format_fun(gwas_non_z)


### PLOT
pA <- plot_fun(ddc_dt, "DDC")
pA
pB <- plot_fun(non_dt, "nonDDC")
pB
pout <- grid.arrange(pA, pB, ncol = 2)


### EXPORT
ggsave("~/Dropbox/PhD/Data/Figures/Figure_SNV_Z_tmp4.jpg", plot = pout, height = 6, width = 12)


