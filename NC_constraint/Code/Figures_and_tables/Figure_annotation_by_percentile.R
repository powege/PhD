rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)

### FUNCTIONS

gobbler <- function(promoter){
  percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                        "51-75", "76-90", "91-95", "96-98", "99", "100")
  mean_fold_change <- rep(NA, 12)
  mean_fold_change[1] <- mean(promoter$fold_change[promoter$percentile %in% c(1)])
  mean_fold_change[2] <- mean(promoter$fold_change[promoter$percentile %in% c(2)])
  mean_fold_change[3] <- mean(promoter$fold_change[promoter$percentile %in% c(3:5)])
  mean_fold_change[4] <- mean(promoter$fold_change[promoter$percentile %in% c(6:10)])
  mean_fold_change[5] <- mean(promoter$fold_change[promoter$percentile %in% c(11:25)])
  mean_fold_change[6] <- mean(promoter$fold_change[promoter$percentile %in% c(26:50)])
  mean_fold_change[7] <- mean(promoter$fold_change[promoter$percentile %in% c(51:75)])
  mean_fold_change[8] <- mean(promoter$fold_change[promoter$percentile %in% c(76:90)])
  mean_fold_change[9] <- mean(promoter$fold_change[promoter$percentile %in% c(91:95)])
  mean_fold_change[10] <- mean(promoter$fold_change[promoter$percentile %in% c(96:98)])
  mean_fold_change[11] <- mean(promoter$fold_change[promoter$percentile %in% c(99)])
  mean_fold_change[12] <- mean(promoter$fold_change[promoter$percentile %in% c(100)])
  out <- data.frame(percentile_group = percentile_group,
                    rank = c(1:12),
                    mean_fold_change = mean_fold_change,
                    annotation = rep(promoter$annotation[1], 12))
  return(out)
}

# FUNCTION to calculate fold change vs the 100th percentile
fold.change <- function(sub){
  sub$fold_change <- sub$fraction/sub$fraction[sub$percentile == 100]
  return(sub)
}

### IMPORT
h_dt_common <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/human_ann_by_percentile_MAF001.csv")
h_dt_all <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/human_ann_by_percentile.csv")
m_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/mouse_ann_by_percentile.csv")

### FORMAT

# subset chr totals
h_dt_all <- subset(h_dt_all, h_dt_all$chromosome == "total")
h_dt_common <- subset(h_dt_common, h_dt_common$chromosome == "total")
m_dt <- subset(m_dt, m_dt$chromosome == "total")
# calculate fraction
h_dt_all$fraction <- h_dt_all$total_annotation_POS/h_dt_all$total_percentile_POS
h_dt_common$fraction <- h_dt_common$total_annotation_POS/h_dt_common$total_percentile_POS
m_dt$fraction <- m_dt$total_annotation_POS/m_dt$total_percentile_POS
# calculate fold change vs the 100th percentile
h_dt_all <- ddply(h_dt_all, "annotation", fold.change)
h_dt_common <- ddply(h_dt_common, "annotation", fold.change)
m_dt <- ddply(m_dt, "annotation", fold.change)
# group percentiles with custom function
h_dt_all <- ddply(h_dt_all, "annotation", gobbler)
h_dt_common <- ddply(h_dt_common, "annotation", gobbler)
m_dt <- ddply(m_dt, "annotation", gobbler)
# set levels order
levels(h_dt_all$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                                   "51-75", "76-90", "91-95", "96-98", "99", "100")
levels(h_dt_common$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                                   "51-75", "76-90", "91-95", "96-98", "99", "100")
levels(m_dt$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                                 "51-75", "76-90", "91-95", "96-98", "99", "100")

### PLOT INDIVIDUAL

m_p <- ggplot(m_dt, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.25, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
m_p
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_mouse_ann_by_percentile.jpg", plot = m_p, height = 5, width = 7)

h_p_all <- ggplot(h_dt_all, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.25, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
h_p_all
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_human_ann_by_percentile_all.jpg", plot = h_p_all, height = 5, width = 7)

h_p_common <- ggplot(h_dt_common, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
  scale_y_continuous(
                      trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.25, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
h_p_common
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_human_ann_by_percentile_MAF001.jpg", plot = h_p_common, height = 5, width = 7)


### PLOT COMBINED

h_dt_common$percentile_group <- factor(h_dt_common$percentile_group,
                                       levels = c("100", "99", "96-98", "91-95", "76-90", "51-75",
                                          "26-50", "11-25", "6-10", "3-5", "2", "1"))
hp1 <- ggplot(h_dt_common, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
  ggtitle("Human") +
  scale_y_continuous( position = "right",
                      trans = log2_trans(),
                      # breaks = trans_breaks("log2", function(x) 2^x)) +
                      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
                      limits = c(0.25, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     trans = "reverse",
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
hp1

mp1 <- ggplot(m_dt, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
  ggtitle("Mouse") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.25, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
mp1

pout <- grid.arrange(mp1, hp1, nrow = 1, widths = c(5, 6.5))
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_HM_ann_by_percentile.jpg", plot = pout, height = 6, width = 10)

# ### FUNCTIONS
# 
# gobbler <- function(promoter){
#   percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                         "51-75", "76-90", "91-95", "96-98", "99", "100")
#   mean_fold_change <- rep(NA, 12)
#   mean_fold_change[1] <- mean(promoter$fold_change[promoter$percentile %in% c(1)])
#   mean_fold_change[2] <- mean(promoter$fold_change[promoter$percentile %in% c(2)])
#   mean_fold_change[3] <- mean(promoter$fold_change[promoter$percentile %in% c(3:5)])
#   mean_fold_change[4] <- mean(promoter$fold_change[promoter$percentile %in% c(6:10)])
#   mean_fold_change[5] <- mean(promoter$fold_change[promoter$percentile %in% c(11:25)])
#   mean_fold_change[6] <- mean(promoter$fold_change[promoter$percentile %in% c(26:50)])
#   mean_fold_change[7] <- mean(promoter$fold_change[promoter$percentile %in% c(51:75)])
#   mean_fold_change[8] <- mean(promoter$fold_change[promoter$percentile %in% c(76:90)])
#   mean_fold_change[9] <- mean(promoter$fold_change[promoter$percentile %in% c(91:95)])
#   mean_fold_change[10] <- mean(promoter$fold_change[promoter$percentile %in% c(96:98)])
#   mean_fold_change[11] <- mean(promoter$fold_change[promoter$percentile %in% c(99)])
#   mean_fold_change[12] <- mean(promoter$fold_change[promoter$percentile %in% c(100)])
#   out <- data.frame(percentile_group = percentile_group,
#                     rank = c(1:12),
#                     mean_fold_change = mean_fold_change,
#                     annotation = rep(promoter$annotation[1], 12))
#   return(out)
# }
# 
# 
# ### IMPORT 
# h_dt_common <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/human_ann_by_percentile_MAF001.csv")
# m_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/mouse_ann_by_percentile.csv")
# 
# ### FORMAT
# 
# # subset chr totals
# h_dt_common <- subset(h_dt_common, h_dt_common$chromosome == "total")
# m_dt <- subset(m_dt, m_dt$chromosome == "total")
# # calculate fraction 
# h_dt_common$fraction <- h_dt_common$total_annotation_POS/h_dt_common$total_percentile_POS
# m_dt$fraction <- m_dt$total_annotation_POS/m_dt$total_percentile_POS
# # calculate fold change vs average percentile
# h_dt_common$fold_change[h_dt_common$annotation == "Enhancer"] <- h_dt_common$fraction[h_dt_common$annotation == "Enhancer"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Enhancer"])
# h_dt_common$fold_change[h_dt_common$annotation == "Exon - CDS"] <- h_dt_common$fraction[h_dt_common$annotation == "Exon - CDS"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Exon - CDS"])
# h_dt_common$fold_change[h_dt_common$annotation == "Exon - non-coding"] <- h_dt_common$fraction[h_dt_common$annotation == "Exon - non-coding"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Exon - non-coding"])
# h_dt_common$fold_change[h_dt_common$annotation == "Exon - UTR"] <- h_dt_common$fraction[h_dt_common$annotation == "Exon - UTR"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Exon - UTR"])
# h_dt_common$fold_change[h_dt_common$annotation == "Intron"] <- h_dt_common$fraction[h_dt_common$annotation == "Intron"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Intron"])
# h_dt_common$fold_change[h_dt_common$annotation == "Open chromatin"] <- h_dt_common$fraction[h_dt_common$annotation == "Open chromatin"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Open chromatin"])
# h_dt_common$fold_change[h_dt_common$annotation == "Promoter"] <- h_dt_common$fraction[h_dt_common$annotation == "Promoter"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Promoter"])
# h_dt_common$fold_change[h_dt_common$annotation == "Promoter flanking"] <- h_dt_common$fraction[h_dt_common$annotation == "Promoter flanking"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Promoter flanking"])
# h_dt_common$fold_change[h_dt_common$annotation == "TF binding"] <- h_dt_common$fraction[h_dt_common$annotation == "TF binding"] / mean(h_dt_common$fraction[h_dt_common$annotation == "TF binding"])
# h_dt_common$fold_change[h_dt_common$annotation == "Unannotated"] <- h_dt_common$fraction[h_dt_common$annotation == "Unannotated"] / mean(h_dt_common$fraction[h_dt_common$annotation == "Unannotated"])
# m_dt$fold_change[m_dt$annotation == "Enhancer"] <- m_dt$fraction[m_dt$annotation == "Enhancer"] / mean(m_dt$fraction[m_dt$annotation == "Enhancer"])
# m_dt$fold_change[m_dt$annotation == "Exon - CDS"] <- m_dt$fraction[m_dt$annotation == "Exon - CDS"] / mean(m_dt$fraction[m_dt$annotation == "Exon - CDS"])
# m_dt$fold_change[m_dt$annotation == "Exon - non-coding"] <- m_dt$fraction[m_dt$annotation == "Exon - non-coding"] / mean(m_dt$fraction[m_dt$annotation == "Exon - non-coding"])
# m_dt$fold_change[m_dt$annotation == "Exon - UTR"] <- m_dt$fraction[m_dt$annotation == "Exon - UTR"] / mean(m_dt$fraction[m_dt$annotation == "Exon - UTR"])
# m_dt$fold_change[m_dt$annotation == "Intron"] <- m_dt$fraction[m_dt$annotation == "Intron"] / mean(m_dt$fraction[m_dt$annotation == "Intron"])
# m_dt$fold_change[m_dt$annotation == "Open chromatin"] <- m_dt$fraction[m_dt$annotation == "Open chromatin"] / mean(m_dt$fraction[m_dt$annotation == "Open chromatin"])
# m_dt$fold_change[m_dt$annotation == "Promoter"] <- m_dt$fraction[m_dt$annotation == "Promoter"] / mean(m_dt$fraction[m_dt$annotation == "Promoter"])
# m_dt$fold_change[m_dt$annotation == "Promoter flanking"] <- m_dt$fraction[m_dt$annotation == "Promoter flanking"] / mean(m_dt$fraction[m_dt$annotation == "Promoter flanking"])
# m_dt$fold_change[m_dt$annotation == "TF binding"] <- m_dt$fraction[m_dt$annotation == "TF binding"] / mean(m_dt$fraction[m_dt$annotation == "TF binding"])
# m_dt$fold_change[m_dt$annotation == "Unannotated"] <- m_dt$fraction[m_dt$annotation == "Unannotated"] / mean(m_dt$fraction[m_dt$annotation == "Unannotated"])
# # # calculate fold change vs the 100th percentile
# # h_dt_common <- ddply(h_dt_common, "annotation", fold.change)
# # m_dt <- ddply(m_dt, "annotation", fold.change)
# # group percentiles with custom function
# h_dt_common <- ddply(h_dt_common, "annotation", gobbler)
# m_dt <- ddply(m_dt, "annotation", gobbler)
# # set levels order
# levels(h_dt_common$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                           "51-75", "76-90", "91-95", "96-98", "99", "100")
# levels(m_dt$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                    "51-75", "76-90", "91-95", "96-98", "99", "100")
# 
# ### PLOT INDIVIDUAL
# 
# m_p <- ggplot(m_dt, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   scale_y_continuous(trans = log2_trans(),
#                      # breaks = trans_breaks("log2", function(x) 2^x)) +
#                      breaks = c(0.15, 0.3, 0.5, 1, 2, 5, 10),
#                      limits = c(0.1, 13)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# m_p
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_mouse_ann_by_percentile.jpg", plot = m_p, height = 5, width = 7)
# 
# h_p_common <- ggplot(h_dt_common, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   scale_y_continuous(
#     trans = log2_trans(),
#     # breaks = trans_breaks("log2", function(x) 2^x)) +
#     breaks = c(0.15, 0.3, 0.5, 1, 2, 5, 12),
#     limits = c(0.1, 13)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# h_p_common
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_human_ann_by_percentile_MAF001.jpg", plot = h_p_common, height = 5, width = 7)
# 
# 
# ### PLOT COMBINED
# 
# h_dt_common$percentile_group <- factor(h_dt_common$percentile_group, 
#                                        levels = c("100", "99", "96-98", "91-95", "76-90", "51-75",
#                                                   "26-50", "11-25", "6-10", "3-5", "2", "1"))
# hp1 <- ggplot(h_dt_common, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   ggtitle("Human") +
#   scale_y_continuous( position = "right",
#                       trans = log2_trans(),
#                       # breaks = trans_breaks("log2", function(x) 2^x)) +
#                       breaks = c(0.15, 0.3, 0.5, 1, 2, 5, 10),
#                       limits = c(0.1, 13)) +
#   scale_x_continuous(breaks=c(1:12),
#                      trans = "reverse",
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# hp1
# 
# mp1 <- ggplot(m_dt, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   ggtitle("Mouse") +
#   scale_y_continuous(trans = log2_trans(),
#                      # breaks = trans_breaks("log2", function(x) 2^x)) +
#                      breaks = c(0.15, 0.3, 0.5, 1, 2, 5, 10),
#                      limits = c(0.1, 13)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# mp1
# 
# pout <- grid.arrange(mp1, hp1, nrow = 1, widths = c(5, 6.5))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_HM_ann_by_percentile.jpg", plot = pout, height = 6, width = 10)



#####

# ### FUNCTIONS
# 
# gobbler <- function(promoter){
#   percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                         "51-75", "76-90", "91-95", "96-98", "99", "100")
#   mean_fold_change <- rep(NA, 12)
#   mean_fold_change[1] <- mean(promoter$fold_change[promoter$percentile %in% c(1)])
#   mean_fold_change[2] <- mean(promoter$fold_change[promoter$percentile %in% c(2)])
#   mean_fold_change[3] <- mean(promoter$fold_change[promoter$percentile %in% c(3:5)])
#   mean_fold_change[4] <- mean(promoter$fold_change[promoter$percentile %in% c(6:10)])
#   mean_fold_change[5] <- mean(promoter$fold_change[promoter$percentile %in% c(11:25)])
#   mean_fold_change[6] <- mean(promoter$fold_change[promoter$percentile %in% c(26:50)])
#   mean_fold_change[7] <- mean(promoter$fold_change[promoter$percentile %in% c(51:75)])
#   mean_fold_change[8] <- mean(promoter$fold_change[promoter$percentile %in% c(76:90)])
#   mean_fold_change[9] <- mean(promoter$fold_change[promoter$percentile %in% c(91:95)])
#   mean_fold_change[10] <- mean(promoter$fold_change[promoter$percentile %in% c(96:98)])
#   mean_fold_change[11] <- mean(promoter$fold_change[promoter$percentile %in% c(99)])
#   mean_fold_change[12] <- mean(promoter$fold_change[promoter$percentile %in% c(100)])
#   out <- data.frame(percentile_group = percentile_group,
#                     rank = c(1:12),
#                     mean_fold_change = mean_fold_change,
#                     annotation = rep(promoter$annotation[1], 12))
#   return(out)
# }
# 
# # FUNCTION to calculate fold change vs the 100th percentile
# fold.change <- function(sub){
#   sub$fold_change <- sub$fraction/sub$fraction[sub$percentile == 100]
#   return(sub)
# }
# 
# ### IMPORT 
# h_dt_common <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/human_ann_by_percentile_MAF001.csv")
# h_dt_all <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/human_ann_by_percentile.csv")
# m_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/mouse_ann_by_percentile.csv")
# 
# ### FORMAT
# 
# # subset chr totals
# h_dt_all <- subset(h_dt_all, h_dt_all$chromosome == "total")
# h_dt_common <- subset(h_dt_common, h_dt_common$chromosome == "total")
# m_dt <- subset(m_dt, m_dt$chromosome == "total")
# # calculate fraction 
# h_dt_all$fraction <- h_dt_all$total_annotation_POS/h_dt_all$total_percentile_POS
# h_dt_common$fraction <- h_dt_common$total_annotation_POS/h_dt_common$total_percentile_POS
# m_dt$fraction <- m_dt$total_annotation_POS/m_dt$total_percentile_POS
# # calculate fold change vs the 100th percentile
# h_dt_all <- ddply(h_dt_all, "annotation", fold.change)
# h_dt_common <- ddply(h_dt_common, "annotation", fold.change)
# m_dt <- ddply(m_dt, "annotation", fold.change)
# # group percentiles with custom function
# h_dt_all <- ddply(h_dt_all, "annotation", gobbler)
# h_dt_common <- ddply(h_dt_common, "annotation", gobbler)
# m_dt <- ddply(m_dt, "annotation", gobbler)
# # set levels order
# levels(h_dt_all$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                    "51-75", "76-90", "91-95", "96-98", "99", "100")
# levels(h_dt_common$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                    "51-75", "76-90", "91-95", "96-98", "99", "100")
# levels(m_dt$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                  "51-75", "76-90", "91-95", "96-98", "99", "100")
# 
# ### PLOT INDIVIDUAL
# 
# m_p <- ggplot(m_dt, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   scale_y_continuous(trans = log2_trans(),
#                      # breaks = trans_breaks("log2", function(x) 2^x)) +
#                      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
#                      limits = c(0.25, 40)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# m_p
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_mouse_ann_by_percentile.jpg", plot = m_p, height = 5, width = 7)
# 
# h_p_all <- ggplot(h_dt_all, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   scale_y_continuous(trans = log2_trans(),
#                      # breaks = trans_breaks("log2", function(x) 2^x)) +
#                      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
#                      limits = c(0.25, 40)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# h_p_all
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_human_ann_by_percentile_all.jpg", plot = h_p_all, height = 5, width = 7)
# 
# h_p_common <- ggplot(h_dt_common, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   scale_y_continuous(
#                       trans = log2_trans(),
#                      # breaks = trans_breaks("log2", function(x) 2^x)) +
#                      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
#                      limits = c(0.25, 40)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# h_p_common
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_human_ann_by_percentile_MAF001.jpg", plot = h_p_common, height = 5, width = 7)
# 
# 
# ### PLOT COMBINED
# 
# h_dt_common$percentile_group <- factor(h_dt_common$percentile_group, 
#                                        levels = c("100", "99", "96-98", "91-95", "76-90", "51-75",
#                                           "26-50", "11-25", "6-10", "3-5", "2", "1"))
# hp1 <- ggplot(h_dt_common, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   ggtitle("Human") +
#   scale_y_continuous( position = "right",
#                       trans = log2_trans(),
#                       # breaks = trans_breaks("log2", function(x) 2^x)) +
#                       breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
#                       limits = c(0.25, 40)) +
#   scale_x_continuous(breaks=c(1:12),
#                      trans = "reverse",
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# hp1
# 
# mp1 <- ggplot(m_dt, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   xlab("Constraint (percentile rank)") +
#   ylab("Genomic territory\n(fold change versus 100th percentile)") +
#   ggtitle("Mouse") +
#   scale_y_continuous(trans = log2_trans(),
#                      # breaks = trans_breaks("log2", function(x) 2^x)) +
#                      breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
#                      limits = c(0.25, 40)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# mp1
# 
# pout <- grid.arrange(mp1, hp1, nrow = 1, widths = c(5, 6.5))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_HM_ann_by_percentile.jpg", plot = pout, height = 6, width = 10)

