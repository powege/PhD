rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(gridExtra)

### PLOT CHANGE RELATIVE TO ANNOTATION AVERAGE 
### RELATIVE!!!!!
### standardise relative to differnce between median and 0

### IMPORT

h_dt_cg <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_CpG_proportion_human.csv")
m_dt_cg <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_CpG_proportion_mouse.csv")

### FORMAT

# add species
h_dt_cg$Species <- "Human"
m_dt_cg$Species <- "Mouse"

# calculate average for all ann
h_ann_cg_median <- median(h_dt_cg$CpG_proportion)
m_ann_cg_median <- median(m_dt_cg$CpG_proportion)
h_ann_cg_mean <- mean(h_dt_cg$CpG_proportion)
m_ann_cg_mean <- mean(m_dt_cg$CpG_proportion)

# calculate deviation from average
# h_dt_cg$CpG_proportion_dif <- (h_dt_cg$CpG_proportion + 1) / (h_ann_cg_median + 1)
# m_dt_cg$CpG_proportion_dif <- (m_dt_cg$CpG_proportion + 1) / (m_ann_cg_median + 1)
h_dt_cg$CpG_proportion_dif <- h_dt_cg$CpG_proportion - h_ann_cg_median
m_dt_cg$CpG_proportion_dif <- m_dt_cg$CpG_proportion - m_ann_cg_median

# rbind
dt_cg <- rbind(h_dt_cg, m_dt_cg)
rm(h_dt_cg, m_dt_cg)

# QC
plot_dt_cg <- dt_cg
# plot_dt_cg <- subset(dt_cg, dt_cg$N_proportion == 0)
# plot_dt_cg <- subset(dt_cg, dt_cg$N_proportion == 0 & dt_cg$length >= 4)

# calculate 1st 2nd and 3rd quartiles for each annotation
tmp_cg1 <- aggregate(plot_dt_cg[, 3], list(plot_dt_cg$Species, plot_dt_cg$category), median)
colnames(tmp_cg1) <- c("Species", "Category", "Q2")
tmp_cg2 <- aggregate(plot_dt_cg[, 3], list(plot_dt_cg$Species, plot_dt_cg$category), quantile, probs = 0.25)
colnames(tmp_cg2) <- c("Species", "Category", "Q1")
tmp_cg3 <- aggregate(plot_dt_cg[, 3], list(plot_dt_cg$Species, plot_dt_cg$category), quantile, probs = 0.75)
colnames(tmp_cg3) <- c("Species", "Category", "Q3")
plot_tmp_cg <- merge(tmp_cg2, tmp_cg1)
plot_tmp_cg <- merge(plot_tmp_cg, tmp_cg3)

# calculate fold change B/A - 1
test <- plot_tmp_cg
test$Q1[test$Species == "Mouse"] <- (test$Q1[test$Species == "Mouse"]/m_ann_cg_median) 
test$Q2[test$Species == "Mouse"] <- (test$Q2[test$Species == "Mouse"]/m_ann_cg_median) 
test$Q3[test$Species == "Mouse"] <- (test$Q3[test$Species == "Mouse"]/m_ann_cg_median) 
test$Q1[test$Species == "Human"] <- (test$Q1[test$Species == "Human"]/h_ann_cg_median) 
test$Q2[test$Species == "Human"] <- (test$Q2[test$Species == "Human"]/h_ann_cg_median) 
test$Q3[test$Species == "Human"] <- (test$Q3[test$Species == "Human"]/h_ann_cg_median) 

# test$Q1[test$Species == "Mouse"] <- (test$Q1[test$Species == "Mouse"] - m_ann_cg_median) / m_ann_cg_median
# test$Q2[test$Species == "Mouse"] <- (test$Q2[test$Species == "Mouse"] - m_ann_cg_median) / m_ann_cg_median
# test$Q3[test$Species == "Mouse"] <- (test$Q3[test$Species == "Mouse"] - m_ann_cg_median) / m_ann_cg_median
# test$Q1[test$Species == "Human"] <- (test$Q1[test$Species == "Human"] - h_ann_cg_median) / m_ann_cg_median
# test$Q2[test$Species == "Human"] <- (test$Q2[test$Species == "Human"] - h_ann_cg_median) / m_ann_cg_median
# test$Q3[test$Species == "Human"] <- (test$Q3[test$Species == "Human"] - h_ann_cg_median) / m_ann_cg_median
plot_tmp_cg <- test

# set factor order by median CpG fraction
order_cg <- aggregate(plot_tmp_cg$Q2~plot_tmp_cg$Category, FUN=median)
order_cg <- order_cg[order(order_cg$`plot_tmp_cg$Q2`),]
plot_tmp_cg$Category <- factor(plot_tmp_cg$Category, levels = as.character(order_cg$`plot_tmp_cg$Category`))
# plot_tmp_cg$Category <- factor(plot_tmp_cg$Category, levels = as.character(order_snv$`plot_tmp_snv$Category`))
# plot_tmp_cg$Species <- factor(plot_tmp_cg$Species, levels = c("Human", "Mouse"))


# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    Category = levels(plot_tmp_cg$Category),
                    rect_type = c("a", "c"))
rects <- rbind(rects)
plot_tmp_cg <- merge(plot_tmp_cg, rects)


### PLOT

p1 <- ggplot(plot_tmp_cg, aes(x=Category, y=Q2, color=Species)) + 
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.5,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  geom_rect(
    aes(xmin = plot_tmp_cg$xmin,
        xmax = plot_tmp_cg$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = plot_tmp_cg$rect_type),
    color = NA,
    alpha = 0.5,
    show.legend = F) +
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  geom_hline(yintercept=1,
             linetype="dashed",
             color = "black",
             size=1) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Genomic annotation") +
  ylab("Fold change") +
  ggtitle("CpG fraction") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8),
                     limits = c(0, 8)) +
  coord_flip() +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    legend.justification=c(1,0),
    legend.position=c(0.95, 0.05),
    legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
p1

### IMPORT

h_dt_snv <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_SNV_proportion_human.csv")
m_dt_snv <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_SNV_proportion_mouse.csv")

### Format

# add species
h_dt_snv$Species <- "Human"
m_dt_snv$Species <- "Mouse"

# calculate average for all ann
h_ann_snv_median <- median(h_dt_snv$n_SNV_kb)
m_ann_snv_median <- median(m_dt_snv$n_SNV_kb)
h_ann_snv_mean <- mean(h_dt_snv$n_SNV_kb)
m_ann_snv_mean <- mean(m_dt_snv$n_SNV_kb)

# calculate deviation from average
h_dt_snv$n_SNV_kb_dif <- h_dt_snv$n_SNV_kb - h_ann_snv_median
m_dt_snv$n_SNV_kb_dif <- m_dt_snv$n_SNV_kb - m_ann_snv_median

# rbind
dt_snv <- rbind(h_dt_snv, m_dt_snv)
rm(h_dt_snv, m_dt_snv)

# QC
plot_dt_snv <- dt_snv
# plot_dt_snv <- subset(dt_snv, dt_snv$N_proportion == 0 & dt_snv$length >= 4)

# calculate 1st 2nd and 3rd quartiles for each annotation
tmp_snv1 <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), median)
colnames(tmp_snv1) <- c("Species", "Category", "Q2")
tmp_snv2 <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), quantile, probs = 0.25)
colnames(tmp_snv2) <- c("Species", "Category", "Q1")
tmp_snv3 <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), quantile, probs = 0.75)
colnames(tmp_snv3) <- c("Species", "Category", "Q3")
plot_tmp_snv <- merge(tmp_snv2, tmp_snv1)
plot_tmp_snv <- merge(plot_tmp_snv, tmp_snv3)

# calculate fold change B/A - 1
test <- plot_tmp_snv
test$Q1[test$Species == "Mouse"] <- (test$Q1[test$Species == "Mouse"]/m_ann_snv_median) 
test$Q2[test$Species == "Mouse"] <- (test$Q2[test$Species == "Mouse"]/m_ann_snv_median) 
test$Q3[test$Species == "Mouse"] <- (test$Q3[test$Species == "Mouse"]/m_ann_snv_median) 
test$Q1[test$Species == "Human"] <- (test$Q1[test$Species == "Human"]/h_ann_snv_median) 
test$Q2[test$Species == "Human"] <- (test$Q2[test$Species == "Human"]/h_ann_snv_median) 
test$Q3[test$Species == "Human"] <- (test$Q3[test$Species == "Human"]/h_ann_snv_median) 

# test$Q1[test$Species == "Mouse"] <- (test$Q1[test$Species == "Mouse"] - m_ann_snv_median) / m_ann_snv_median
# test$Q2[test$Species == "Mouse"] <- (test$Q2[test$Species == "Mouse"] - m_ann_snv_median) / m_ann_snv_median
# test$Q3[test$Species == "Mouse"] <- (test$Q3[test$Species == "Mouse"] - m_ann_snv_median) / m_ann_snv_median
# test$Q1[test$Species == "Human"] <- (test$Q1[test$Species == "Human"] - h_ann_snv_median) / m_ann_snv_median
# test$Q2[test$Species == "Human"] <- (test$Q2[test$Species == "Human"] - h_ann_snv_median) / m_ann_snv_median
# test$Q3[test$Species == "Human"] <- (test$Q3[test$Species == "Human"] - h_ann_snv_median) / m_ann_snv_median
plot_tmp_snv <- test

# set factor order by median CpG fraction
order_snv <- aggregate(plot_tmp_snv$Q2~plot_tmp_snv$Category, FUN=median)
order_snv <- order_snv[order(order_snv$`plot_tmp_snv$Q2`),]
plot_tmp_snv$Category <- factor(plot_tmp_snv$Category, levels = as.character(order_cg$`plot_tmp_cg$Category`))
# plot_tmp_snv$Category <- factor(plot_tmp_snv$Category, levels = as.character(order_snv$`plot_tmp_snv$Category`))
# plot_tmp_snv$Species <- factor(plot_tmp_snv$Species, levels = c("Human", "Mouse"))


# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    Category = levels(plot_tmp_snv$Category),
                    rect_type = c("a", "c"))
rects <- rbind(rects)
plot_tmp_snv <- merge(plot_tmp_snv, rects)


### PLOT

p2 <- ggplot(plot_tmp_snv, aes(x=Category, y=Q2, color=Species)) + 
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.5,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  geom_rect(
    aes(xmin = plot_tmp_snv$xmin,
        xmax = plot_tmp_snv$xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = plot_tmp_snv$rect_type),
    color = NA,
    alpha = 0.5,
    show.legend = F) +
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  geom_hline(yintercept=1, 
             linetype="dashed", 
             color = "black", 
             size=1) +
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Genomic annotation") +
  ylab("Fold change") +
  ggtitle("N SNVs per kb") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2),
                     limits = c(0, 2)) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.title = element_blank(),
    # legend.justification=c(1,0),
    # legend.position=c(0.23, 0.70),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
p2

pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 6))
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ann_CpG_SNV_prop.jpg", plot = pout, height = 6, width = 10)


# ### IMPORT
# 
# h_dt_cg <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_CpG_proportion_human.csv")
# m_dt_cg <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_CpG_proportion_mouse.csv")
# 
# ### FORMAT
# 
# # add species
# h_dt_cg$Species <- "Human"
# m_dt_cg$Species <- "Mouse"
# 
# # calculate average for all ann
# h_ann_cg_median <- median(h_dt_cg$CpG_proportion)
# m_ann_cg_median <- median(m_dt_cg$CpG_proportion)
# h_ann_cg_mean <- mean(h_dt_cg$CpG_proportion)
# m_ann_cg_mean <- mean(m_dt_cg$CpG_proportion)
# 
# # rbind
# dt_cg <- rbind(h_dt_cg, m_dt_cg)
# rm(h_dt_cg, m_dt_cg)
# 
# # QC
# # plot_dt_cg <- dt_cg
# # plot_dt_cg <- subset(dt_cg, dt_cg$N_proportion == 0)
# # plot_dt_cg <- subset(dt_cg, dt_cg$N_proportion == 0 & dt_cg$length >= 4)
# 
# # calculate 1st 2nd and 3rd quartiles for each annotation
# tmp_cg1 <- aggregate(plot_dt_cg[, 3], list(plot_dt_cg$Species, plot_dt_cg$category), median)
# colnames(tmp_cg1) <- c("Species", "Category", "Q2")
# tmp_cg2 <- aggregate(plot_dt_cg[, 3], list(plot_dt_cg$Species, plot_dt_cg$category), quantile, probs = 0.25)
# colnames(tmp_cg2) <- c("Species", "Category", "Q1")
# tmp_cg3 <- aggregate(plot_dt_cg[, 3], list(plot_dt_cg$Species, plot_dt_cg$category), quantile, probs = 0.75)
# colnames(tmp_cg3) <- c("Species", "Category", "Q3")
# plot_tmp_cg <- merge(tmp_cg2, tmp_cg1)
# plot_tmp_cg <- merge(plot_tmp_cg, tmp_cg3)
# 
# # set factor order by median CpG fraction
# order_cg <- aggregate(plot_tmp_cg$Q2~plot_tmp_cg$Category, FUN=median)
# order_cg <- order_cg[order(order_cg$`plot_tmp_cg$Q2`),]
# plot_tmp_cg$Category <- factor(plot_tmp_cg$Category, levels = as.character(order_cg$`plot_tmp_cg$Category`))
# # plot_tmp_cg$Category <- factor(plot_tmp_cg$Category, levels = as.character(order_snv$`plot_tmp_snv$Category`))
# # plot_tmp_cg$Species <- factor(plot_tmp_cg$Species, levels = c("Human", "Mouse"))
# 
# 
# # set rectangle coordinates
# rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1),
#                     xmax = tail(seq, -1),
#                     Category = levels(plot_tmp_cg$Category),
#                     rect_type = c("a", "c"))
# rects <- rbind(rects)
# plot_tmp_cg <- merge(plot_tmp_cg, rects)
# 
# 
# ### PLOT
# 
# p1 <- ggplot(plot_tmp_cg, aes(x=Category, y=Q2, color=Species)) +
#   geom_errorbar(aes(ymax = Q3, ymin = Q1),
#                 position = position_dodge(width=0.9),
#                 stat = "identity",
#                 width=0.5,
#                 size=1.2) +
#   geom_point(position = position_dodge(0.9), size = 2.5) +
#   geom_rect(
#     aes(xmin = plot_tmp_cg$xmin,
#         xmax = plot_tmp_cg$xmax,
#         ymin = -Inf,
#         ymax = Inf,
#         fill = plot_tmp_cg$rect_type),
#     color = NA,
#     alpha = 0.5,
#     show.legend = F) +
#   geom_errorbar(aes(ymax = Q3, ymin = Q1),
#                 position = position_dodge(width=0.9),
#                 stat = "identity",
#                 width=0.4,
#                 size=1.2) +
#   geom_point(position = position_dodge(0.9), size = 2.5) +
#   geom_hline(yintercept=h_ann_cg_median,
#              linetype="dashed",
#              color = "red",
#              size=1) +
#   geom_hline(yintercept=m_ann_cg_median,
#              linetype="dashed",
#              color = "blue",
#              size=1) +
#   scale_fill_manual(values = c("grey", "white")) +
#   xlab("Genomic annotation") +
#   ylab("CpG fraction") +
#   scale_y_continuous(breaks = c(0, 0.04, 0.08, 0.12, 0.16),
#                      limits = c(0, 0.165)) +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     # legend.position = "none",
#     legend.title = element_blank(),
#     legend.justification=c(1,0),
#     legend.position=c(0.95, 0.05),
#     legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14)
#   )
# p1
# 
# ### IMPORT
# 
# h_dt_snv <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_SNV_proportion_human.csv")
# m_dt_snv <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_SNV_proportion_mouse.csv")
# 
# ### Format
# 
# # add species
# h_dt_snv$Species <- "Human"
# m_dt_snv$Species <- "Mouse"
# 
# # calculate average for all ann
# h_ann_snv_median <- median(h_dt_snv$n_SNV_kb)
# m_ann_snv_median <- median(m_dt_snv$n_SNV_kb)
# h_ann_snv_mean <- mean(h_dt_snv$n_SNV_kb)
# m_ann_snv_mean <- mean(m_dt_snv$n_SNV_kb)
# 
# # rbind
# dt_snv <- rbind(h_dt_snv, m_dt_snv)
# rm(h_dt_snv, m_dt_snv)
# 
# # QC
# plot_dt_snv <- dt_snv
# # plot_dt_snv <- subset(dt_snv, dt_snv$N_proportion == 0 & dt_snv$length >= 4)
# 
# # calculate 1st 2nd and 3rd quartiles for each annotation
# tmp_snv1 <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), median)
# colnames(tmp_snv1) <- c("Species", "Category", "Q2")
# tmp_snv2 <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), quantile, probs = 0.25)
# colnames(tmp_snv2) <- c("Species", "Category", "Q1")
# tmp_snv3 <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), quantile, probs = 0.75)
# colnames(tmp_snv3) <- c("Species", "Category", "Q3")
# plot_tmp_snv <- merge(tmp_snv2, tmp_snv1)
# plot_tmp_snv <- merge(plot_tmp_snv, tmp_snv3)
# 
# # set factor order by median CpG fraction
# order_snv <- aggregate(plot_tmp_snv$Q2~plot_tmp_snv$Category, FUN=median)
# order_snv <- order_snv[order(order_snv$`plot_tmp_snv$Q2`),]
# plot_tmp_snv$Category <- factor(plot_tmp_snv$Category, levels = as.character(order_cg$`plot_tmp_cg$Category`))
# # plot_tmp_snv$Category <- factor(plot_tmp_snv$Category, levels = as.character(order_snv$`plot_tmp_snv$Category`))
# # plot_tmp_snv$Species <- factor(plot_tmp_snv$Species, levels = c("Human", "Mouse"))
# 
# 
# # set rectangle coordinates
# rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1),
#                     xmax = tail(seq, -1),
#                     Category = levels(plot_tmp_snv$Category),
#                     rect_type = c("a", "c"))
# rects <- rbind(rects)
# plot_tmp_snv <- merge(plot_tmp_snv, rects)
# 
# 
# ### PLOT
# 
# p2 <- ggplot(plot_tmp_snv, aes(x=Category, y=Q2, color=Species)) +
#   geom_errorbar(aes(ymax = Q3, ymin = Q1),
#                 position = position_dodge(width=0.9),
#                 stat = "identity",
#                 width=0.5,
#                 size=1.2) +
#   geom_point(position = position_dodge(0.9), size = 2.5) +
#   geom_rect(
#     aes(xmin = plot_tmp_snv$xmin,
#         xmax = plot_tmp_snv$xmax,
#         ymin = -Inf,
#         ymax = Inf,
#         fill = plot_tmp_snv$rect_type),
#     color = NA,
#     alpha = 0.5,
#     show.legend = F) +
#   geom_errorbar(aes(ymax = Q3, ymin = Q1),
#                 position = position_dodge(width=0.9),
#                 stat = "identity",
#                 width=0.4,
#                 size=1.2) +
#   geom_point(position = position_dodge(0.9), size = 2.5) +
#   geom_hline(yintercept=h_ann_snv_median,
#              linetype="dashed",
#              color = "red",
#              size=1) +
#   geom_hline(yintercept=m_ann_snv_median,
#              linetype="dashed",
#              color = "blue",
#              size=1) +
#   scale_fill_manual(values = c("grey", "white")) +
#   xlab("Genomic annotation") +
#   ylab("N SNVs per kb") +
#   scale_y_continuous(breaks = c(0, 10, 20, 30, 40),
#                      limits = c(0, 40)) +
#   scale_x_discrete(position = "top") +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     # legend.title = element_blank(),
#     # legend.justification=c(1,0),
#     # legend.position=c(0.23, 0.70),
#     # legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14)
#   )
# p2
# 
# pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 6))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ann_CpG_SNV_prop.jpg", plot = pout, height = 6, width = 10)


######

# xM <- subset(dt_snv, dt_snv$Species == "Mouse")
# xM <- subset(xM, xM$n_SNV == 0)
# table(xM$category)
# 
# xH <- subset(dt_snv, dt_snv$Species == "Human")
# xH <- subset(xH, xH$n_SNV == 0)
# table(xH$category)
# 
# summary(dt_snv$length_kb[dt_snv$Species == "Mouse"])
# test <- aggregate(plot_dt_snv[, 5], list(plot_dt_snv$Species, plot_dt_snv$category), median)

