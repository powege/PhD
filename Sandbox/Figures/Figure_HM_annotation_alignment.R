rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)



dt_plot2 <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_alignment_by_human_annotation.csv")
dt_plot <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/HM_alignment_annotation_summary.csv")

# format total alignment
tmp_tot <- dt_plot[,c("H_annotation", "H_ann_align_frac")]
colnames(tmp_tot) <- c("category", "Total_aligned")
tmp_tot <- tmp_tot[!duplicated(tmp_tot),]

# set factor order by total alignment
# order_H <- aggregate(dt_plot2$align_frac~dt_plot2$category, FUN=median)
order_H <- aggregate(dt_plot$H_ann_align_frac~dt_plot$H_annotation, FUN=median)
colnames(order_H) <- c("category", "median_align")
order_H <- order_H[order(order_H$median_align),]
order_H <- order_H$category
order_M <- rev(order_H)
dt_plot2$category <- factor(dt_plot2$category, levels = as.character(order_H))
dt_plot$M_annotation <- factor(dt_plot$M_annotation, levels = as.character(order_M))
dt_plot$H_annotation <- factor(dt_plot$H_annotation, levels = as.character(order_H))

# calculate 1st 2nd and 3rd quartiles for each annotation
tmp_2 <- aggregate(dt_plot2[, "align_frac"], list(dt_plot2$category), median)
colnames(tmp_2) <- c("category", "Q2")
tmp_1 <- aggregate(dt_plot2[, "align_frac"], list(dt_plot2$category), quantile, probs = 0.25)
colnames(tmp_1) <- c("category", "Q1")
tmp_3 <- aggregate(dt_plot2[, "align_frac"], list(dt_plot2$category), quantile, probs = 0.75)
colnames(tmp_3) <- c("category", "Q3")
plot_tmp <- merge(tmp_2, tmp_1)
plot_tmp <- merge(plot_tmp, tmp_3)
plot_tmp <- merge(plot_tmp, tmp_tot)

# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    category = levels(plot_tmp$category),
                    rect_type = c("a", "c"))
rects <- rbind(rects)
plot_tmp <- merge(plot_tmp, rects)

p2 <- ggplot(plot_tmp, aes(x=category, y=Q2)) + 
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(size = 2.5) + 
  geom_rect(
    aes(xmin = plot_tmp$xmin,
        xmax = plot_tmp$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = plot_tmp$rect_type),
    color = NA,
    alpha = 0.5,
    show.legend = F) +
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(size = 2.5) + 
  geom_point(aes(x=category, y=Total_aligned), colour="red", shape=4, size=3, stroke=2) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment (fraction)") +
  ggtitle("A") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(0, 1)) +
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
p2
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp2.jpg", plot = p2, height = 6, width = 6)

# p2 <- ggplot(dt_plot2, aes(x=category, y=align_frac, fill=category)) + 
#   geom_boxplot(outlier.colour = NA) +
#   # geom_point(size = 2.5) + 
#   # geom_point(aes(x=Category, y=Total_synteny), colour="red", shape=4, size=3, stroke=2) +
#   xlab("Human genomic annotation") +
#   ylab("Human-Mouse alignment (fraction)") +
#   # ggtitle("") +
#   scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
#                      limits = c(0, 1)) +
#   scale_fill_brewer(palette="Set3", direction = -1) +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     # legend.title = element_blank(),
#     # legend.justification=c(1,0),
#     # legend.position=c(0.95, 0.05),
#     # legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14)
#   )
# p2
# ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp2.jpg", plot = p2, height = 6, width = 6)

p1 <- ggplot() +
  geom_bar(data=dt_plot, aes(x=H_annotation, y=HM_ann_align_frac, fill=M_annotation), colour = "black", stat="identity") +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment composition (fraction)") +
  ggtitle("B") +
  labs(fill = "Mouse genomic\nannotation") +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8),
                     limits = c(0, 1.8)) +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(
    # legend.position = "none",
    # legend.title = element_blank(),
    legend.title = element_text(hjust = 0.5),
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
p1
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp.jpg", plot = p1, height = 6, width = 8)

pout <- grid.arrange(p2, p1, nrow = 1, widths = c(6, 8))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp3.jpg", plot = pout, height = 6, width = 14)





####

# ddply(dt_plot, "N_M_ann_align_frac", sum)
# 
# out_vec1 <- rep(NA, 10)
# out_vec2 <- rep(NA, 10)
# out_vec3 <- rep(NA, 10)
# for (i in 1:length(unique(dt_plot$H_ann))){
#   sub <- subset(dt_plot, dt_plot$H_ann == unique(dt_plot$H_ann)[i])
#   out_vec1[i] <- sum(sub$N_M_ann_align_frac)
#   out_vec2[i] <- sum(sub$N_M_ann_align)
#   out_vec3[i] <- sub$N_H_ann_align[1]
# }

# ggplot(data=dt_plot, aes(x=H_ann, y=N_M_ann_align_frac, colour=M_ann)) +
#   geom_point(shape=4, size=3, stroke=2) +
#   xlab("Human genomic annotation") +
#   ylab("Fraction of orthologous bases in mouse") +
#   coord_flip() +
#   scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8),
#                      limits = c(0.1, 0.9)) 
# tmp <- data.table(out_vec1, out_vec2, out_vec3)
# tmp$frac <- tmp$out_vec2/tmp$out_vec3


# ### CORRELATION PLOT
# 
# library(corrplot)
# library(tidyr)
# library(RColorBrewer)
# 
# cor_dt <- dt_plot[,c("H_annotation", "M_annotation", "HM_ann_align_frac")]
# cor_dt <- spread(cor_dt, M_annotation, HM_ann_align_frac)
# cor_dt <- cor_dt[order(-H_annotation),]
# cor_dt$H_annotation <- NULL
# cor_dt <- as.matrix(cor_dt)
# rownames(cor_dt) <- colnames(cor_dt)
# pdf(file = "~/Dropbox/PhD/Data/Figures/Figure_HM_alignment_matrix.pdf", width = 6, height = 6)
# corrplot(cor_dt,
#          method = "color",
#          cl.lim = c(0, 1),
#          col = brewer.pal(n = 10, name = "RdYlBu"),
#          number.cex = 0.9,
#          addCoef.col = "black",
#          tl.col = "black")
# dev.off()



