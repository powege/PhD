rm(list = ls())
graphics.off()

library(ggplot2)
library(data.table)
library(scales)

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
                    annotation = rep(h_gerp$annotation[1], 12),
                    grouping = rep(h_gerp$grouping[1], 12))
  return(out)
}

### IMPORT
dt <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/ClinVar_SNVs_by_percentile.csv")

### FORMAT
dt$grouping <- paste0(dt$species, " ", dt$annotation)

# calculate fold change vs the 100th percentile
dt$fold_change <- NA
dt$fold_change[dt$grouping == "Human Non-coding"] <- dt$n_ClinVar_adj[dt$grouping == "Human Non-coding"] / dt$n_ClinVar_adj[dt$grouping == "Human Non-coding" & dt$percentile == 100]
dt$fold_change[dt$grouping == "Human Protein-coding"] <- dt$n_ClinVar_adj[dt$grouping == "Human Protein-coding"] / dt$n_ClinVar_adj[dt$grouping == "Human Protein-coding" & dt$percentile == 100]
dt$fold_change[dt$grouping == "Mouse Non-coding"] <- dt$n_ClinVar_adj[dt$grouping == "Mouse Non-coding"] / dt$n_ClinVar_adj[dt$grouping == "Mouse Non-coding" & dt$percentile == 100]
dt$fold_change[dt$grouping == "Mouse Protein-coding"] <- dt$n_ClinVar_adj[dt$grouping == "Mouse Protein-coding"] / dt$n_ClinVar_adj[dt$grouping == "Mouse Protein-coding" & dt$percentile == 100]

# group percentiles with custom function
dt_plot <- rbind(
  gobbler(dt[dt$grouping == "Human Non-coding",]),
  gobbler(dt[dt$grouping == "Human Protein-coding",]),
  gobbler(dt[dt$grouping == "Mouse Non-coding",]),
  gobbler(dt[dt$grouping == "Mouse Protein-coding",])
)

### PLOT 1

p_1 <- ggplot(dt_plot, aes(x=rank, y=mean_fold_change, color=grouping)) + 
  # geom_point() +
  geom_line(size=1.5) +
  geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change), 
                stat = "identity",
                width=0.2,
                size=1) + 
  xlab("Constraint (percentile rank)") +
  ylab("N pathogenic SNVs\n(fold change versus 100th percentile)") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(1, 2, 5, 10, 20, 40),
                     limits = c(0.5, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.justification=c(1,0), 
    legend.position=c(0.95, 0.70),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_1
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_ClinVar_by_percentile.jpg", plot = p_1, height = 5, width = 6)


### PLOT

p_2 <- ggplot(dt, aes(x=percentile, y=n_ClinVar, color=grouping)) + 
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("N ClinVar SNVs") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p_2
# ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_mouse_ann_by_percentile.jpg", plot = m_p, height = 5, width = 7)

