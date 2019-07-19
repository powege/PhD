rm(list = ls())
graphics.off()

library(ggplot2)
library(data.table)
library(scales)
library(gridExtra)

### FUNCTIONS

gobbler <- function(h_gerp){
  percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                        "51-75", "76-90", "91-95", "96-98", "99", "100")
  
  mean_fold_change <- rep(NA, 12)
  mean_fold_change[1] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(1)])
  mean_fold_change[2] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(2)])
  mean_fold_change[3] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])
  mean_fold_change[4] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])
  mean_fold_change[5] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])
  mean_fold_change[6] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])
  mean_fold_change[7] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])
  mean_fold_change[8] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])
  mean_fold_change[9] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])
  mean_fold_change[10] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])
  mean_fold_change[11] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(99)])
  mean_fold_change[12] <- mean(h_gerp$fold_change[h_gerp$percentile %in% c(100)])
  
  upper_range_fold_change <- rep(NA, 12)
  upper_range_fold_change[1] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(1)])[1]
  upper_range_fold_change[2] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(2)])[1]
  upper_range_fold_change[3] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])[1]
  upper_range_fold_change[4] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])[1]
  upper_range_fold_change[5] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])[1]
  upper_range_fold_change[6] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])[1]
  upper_range_fold_change[7] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])[1]
  upper_range_fold_change[8] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])[1]
  upper_range_fold_change[9] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])[1]
  upper_range_fold_change[10] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])[1]
  upper_range_fold_change[11] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(99)])[1]
  upper_range_fold_change[12] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(100)])[1]
  
  lower_range_fold_change <- rep(NA, 12)
  lower_range_fold_change[1] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(1)])[2]
  lower_range_fold_change[2] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(2)])[2]
  lower_range_fold_change[3] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(3:5)])[2]
  lower_range_fold_change[4] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(6:10)])[2]
  lower_range_fold_change[5] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(11:25)])[2]
  lower_range_fold_change[6] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(26:50)])[2]
  lower_range_fold_change[7] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(51:75)])[2]
  lower_range_fold_change[8] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(76:90)])[2]
  lower_range_fold_change[9] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(91:95)])[2]
  lower_range_fold_change[10] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(96:98)])[2]
  lower_range_fold_change[11] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(99)])[2]
  lower_range_fold_change[12] <- range(h_gerp$fold_change[h_gerp$percentile %in% c(100)])[2]
  
  out <- data.frame(percentile_group = percentile_group,
                    rank = c(1:12),
                    mean_fold_change = mean_fold_change,
                    lower_range_fold_change = lower_range_fold_change,
                    upper_range_fold_change = upper_range_fold_change,
                    species = rep(h_gerp$species[1], 12),
                    annotation = rep(h_gerp$annotation[1], 12))
  return(out)
}

### IMPORT
dt1 <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/ClinVar_SNVs_by_percentile.csv")
dt2 <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/GWAS_SNVs_by_percentile.csv")

### FORMAT

dt1$annotation <- paste0("ClinVar ", dt1$annotation)
dt2$annotation <- "GWAS"
colnames(dt1) <- c("percentile", "n_SNV", "fraction", "species", "annotation", "n_SNV_adj")
colnames(dt2) <- c("percentile", "n_SNV", "fraction", "species", "n_SNV_adj", "annotation")
dt <- rbind(dt1, dt2)
h_dt <- subset(dt, dt$species == "Human")
m_dt <- subset(dt, dt$species == "Mouse")

# calculate fold change vs the 100th percentile
h_dt$fold_change[h_dt$annotation == "ClinVar Non-coding"] <- h_dt$n_SNV_adj[h_dt$annotation == "ClinVar Non-coding"] / h_dt$n_SNV_adj[h_dt$annotation == "ClinVar Non-coding" & h_dt$percentile == 100]
h_dt$fold_change[h_dt$annotation == "ClinVar Protein-coding"] <- h_dt$n_SNV_adj[h_dt$annotation == "ClinVar Protein-coding"] / h_dt$n_SNV_adj[h_dt$annotation == "ClinVar Protein-coding" & h_dt$percentile == 100]
h_dt$fold_change[h_dt$annotation == "GWAS"] <- h_dt$n_SNV_adj[h_dt$annotation == "GWAS"] / h_dt$n_SNV_adj[h_dt$annotation == "GWAS" & h_dt$percentile == 100]
m_dt$fold_change[m_dt$annotation == "ClinVar Non-coding"] <- m_dt$n_SNV_adj[m_dt$annotation == "ClinVar Non-coding"] / m_dt$n_SNV_adj[m_dt$annotation == "ClinVar Non-coding" & m_dt$percentile == 100]
m_dt$fold_change[m_dt$annotation == "ClinVar Protein-coding"] <- m_dt$n_SNV_adj[m_dt$annotation == "ClinVar Protein-coding"] / m_dt$n_SNV_adj[m_dt$annotation == "ClinVar Protein-coding" & m_dt$percentile == 100]
m_dt$fold_change[m_dt$annotation == "GWAS"] <- m_dt$n_SNV_adj[m_dt$annotation == "GWAS"] / m_dt$n_SNV_adj[m_dt$annotation == "GWAS" & m_dt$percentile == 100]

h_dt2 <- h_dt
m_dt2 <- m_dt
h_dt2$fold_change <- h_dt2$n_SNV_adj
m_dt2$fold_change <- m_dt2$n_SNV_adj

# group percentiles with custom function
h_dt_plot <- rbind(
  gobbler(h_dt[h_dt$annotation == "ClinVar Non-coding",]),
  gobbler(h_dt[h_dt$annotation == "ClinVar Protein-coding",]),
  gobbler(h_dt[h_dt$annotation == "GWAS",])
)
m_dt_plot <- rbind(
  gobbler(m_dt[m_dt$annotation == "ClinVar Non-coding",]),
  gobbler(m_dt[m_dt$annotation == "ClinVar Protein-coding",]),
  gobbler(m_dt[m_dt$annotation == "GWAS",])
)
h_dt2_plot <- rbind(
  gobbler(h_dt2[h_dt2$annotation == "ClinVar Non-coding",]),
  gobbler(h_dt2[h_dt2$annotation == "ClinVar Protein-coding",]),
  gobbler(h_dt2[h_dt2$annotation == "GWAS",])
)
m_dt2_plot <- rbind(
  gobbler(m_dt2[m_dt2$annotation == "ClinVar Non-coding",]),
  gobbler(m_dt2[m_dt2$annotation == "ClinVar Protein-coding",]),
  gobbler(m_dt2[m_dt2$annotation == "GWAS",])
)

### PLOT RELATIVE

p_m <- ggplot(m_dt_plot, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change),
                stat = "identity",
                width=0.2,
                size=1) +
  xlab("Constraint (percentile rank)") +
  ylab("N SNVs\n(fold change versus 100th percentile)") +
  ggtitle("Mouse") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.05, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.70),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_m


p_h <- ggplot(h_dt_plot, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change),
                stat = "identity",
                width=0.2,
                size=1) +
  xlab("Constraint (percentile rank)") +
  ylab("N SNVs\n(fold change versus 100th percentile)") +
  ggtitle("Human") +
  scale_y_continuous( position = "right",
                      trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.05, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     trans = "reverse",
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    # axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.justification=c(1,0),
    legend.position=c(0.23, 0.70),
    legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_h

pout <- grid.arrange(p_m, p_h, nrow = 1, widths = c(5, 5))
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_HM_SNV_by_percentile.jpg", plot = pout, height = 6, width = 10)

### PLOT ACTUAL

# p_m2 <- ggplot(m_dt2_plot, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change), 
#                 stat = "identity",
#                 width=0.2,
#                 size=1) + 
#   xlab("Constraint (percentile rank)") +
#   ylab("N SNVs\n(fold change versus 100th percentile)") +
#   ggtitle("Mouse") +
#   scale_y_continuous(
#     # breaks = trans_breaks("log2", function(x) 2^x)) +
#     # breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40),
#     limits = c(0, 3300)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     # legend.justification=c(1,0), 
#     # legend.position=c(0.95, 0.70),  
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# p_m2
# 
# 
# p_h2 <- ggplot(h_dt2_plot, aes(x=rank, y=mean_fold_change, color=annotation)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change), 
#                 stat = "identity",
#                 width=0.2,
#                 size=1) + 
#   xlab("Constraint (percentile rank)") +
#   ylab("N SNVs\n(fold change versus 100th percentile)") +
#   ggtitle("Human") +
#   scale_y_continuous( position = "right",
#                       limits = c(0, 3300)) +
#   scale_x_continuous(breaks=c(1:12),
#                      trans = "reverse",
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     # axis.title.y = element_blank(),
#     legend.title = element_blank(),
#     legend.justification=c(1,0), 
#     legend.position=c(0.23, 0.70),  
#     legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# p_h2
# 
# pout <- grid.arrange(p_m2, p_h2, nrow = 1, widths = c(5, 5))
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_HM_SNV_by_percentile_2.jpg", plot = pout, height = 6, width = 10)

m_dt_plot <- subset(m_dt_plot, m_dt_plot$annotation != "GWAS")
p_m <- ggplot(m_dt_plot, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change),
                stat = "identity",
                width=0.2,
                size=1) +
  xlab("Constraint (percentile rank)") +
  ylab("N SNVs\n(fold change versus 100th percentile)") +
  ggtitle("Mouse") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.2, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.70),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_m

h_dt_plot <- subset(h_dt_plot, h_dt_plot$annotation != "GWAS")
p_h <- ggplot(h_dt_plot, aes(x=rank, y=mean_fold_change, color=annotation)) +
  # geom_point() +
  geom_line(size=1.5) +
  geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change),
                stat = "identity",
                width=0.2,
                size=1) +
  xlab("Constraint (percentile rank)") +
  ylab("N SNVs\n(fold change versus 100th percentile)") +
  ggtitle("Human") +
  scale_y_continuous( position = "right",
                      trans = log2_trans(),
                      # breaks = trans_breaks("log2", function(x) 2^x)) +
                      breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 40),
                      limits = c(0.2, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     trans = "reverse",
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    # axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.justification=c(1,0),
    legend.position=c(0.23, 0.70),
    legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_h

pout <- grid.arrange(p_m, p_h, nrow = 1, widths = c(5, 5))
ggsave("~/Dropbox/Figure_HM_SNV_by_percentile.jpg", plot = pout, height = 6, width = 10)

