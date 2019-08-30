rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

tmp <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation.csv")
tmp_tot <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/Human_mouse_synteny_by_human_annotation_totals.csv")

# calculate 1st 2nd and 3rd quartiles for each annotation
tmp_2 <- aggregate(tmp[, "align_frac"], list(tmp$category), median)
colnames(tmp_2) <- c("Category", "Q2")
tmp_1 <- aggregate(tmp[, "align_frac"], list(tmp$category), quantile, probs = 0.25)
colnames(tmp_1) <- c("Category", "Q1")
tmp_3 <- aggregate(tmp[, "align_frac"], list(tmp$category), quantile, probs = 0.75)
colnames(tmp_3) <- c("Category", "Q3")
plot_tmp <- merge(tmp_2, tmp_1)
plot_tmp <- merge(plot_tmp, tmp_3)

# format total synteny
tmp_tot <- tmp_tot[,c("annotation", "ann_aligned_frac")]
colnames(tmp_tot) <- c("Category", "Total_synteny")
plot_tmp <- merge(plot_tmp, tmp_tot)

# set factor order by total synteny
order_cg <- aggregate(plot_tmp$Q2~plot_tmp$Category, FUN=median)
order_cg <- plot_tmp[order(plot_tmp$Total_synteny),]
plot_tmp$Category <- factor(plot_tmp$Category, levels = as.character(order_cg$Category))

# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    Category = levels(plot_tmp$Category),
                    rect_type = c("a", "c"))
rects <- rbind(rects)
plot_tmp <- merge(plot_tmp, rects)


### PLOT

p1 <- ggplot(plot_tmp, aes(x=Category, y=Q2)) + 
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
  geom_point(aes(x=Category, y=Total_synteny), colour="red", shape=4, size=3, stroke=2) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Human genomic annotation") +
  ylab("Mouse synteny (fraction)") +
  # ggtitle("") +
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
    text = element_text(size=14)
  )
p1
ggsave("~/Dropbox/Figure_synteny_by_annotation.jpg", plot = p1, height = 6, width = 6)

### CLinVar
CV <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_mouse_syntenic.csv")


############


p1 <- ggplot(tmp, aes(x=category, y=align_frac, fill=category)) +
  geom_boxplot() +
  coord_flip()
p1


p2 <- ggplot(tmp, aes(x=category, y=align_frac, fill=category)) +
  geom_violin(trim=T) +
  coord_flip() 
  # stat_summary(fun.data=mean_sdl, 
  #              geom="pointrange", 
  #              color="black")
p2

table(tmp$category)
hist(tmp$align_frac[tmp$category == "Enhancer"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Exon - CDS"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Exon - non-coding"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Exon - UTR"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Intron"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Open chromatin"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Promoter"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Promoter flanking"], breaks = 100)
hist(tmp$align_frac[tmp$category == "TF binding"], breaks = 100)
hist(tmp$align_frac[tmp$category == "Unannotated"], breaks = 100)

tmp$Category <- factor(tmp$Category, levels = as.character(order_cg$Category))
ggplot(data=tmp, aes(x=align_frac, group = category, fill = category)) + 
  geom_density(adjust=1.5, position="fill")

ggplot(data=tmp, aes(x=align_frac, group = category, fill = category)) + 
  geom_density(adjust=1.5 , alpha=0.2)
