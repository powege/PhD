rm(list=ls())
graphics.off()

library(ggplot2)
library(data.table)
library(scales)
library(reshape2)
library(plyr)

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
                    annotation = rep(h_gerp$variable[1], 12))
  return(out)
}

### IMPORT 
h_con <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_conservation_by_percentile_human.csv")
m_con <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_conservation_by_percentile_mouse.csv")


### FORMAT

# subset chr totals
h_con <- subset(h_con, h_con$chromosome == "total")
m_con <- subset(m_con, m_con$chromosome == "total")

# add species
m_con$species <- "Mouse"
h_con$species <- "Human"

# rbind 
con <- rbind(h_con, m_con)

# melt 
con <- melt(data = con, id.vars = c("percentile", "chromosome", "total_percentile_POS", "species"),
             measure.vars = c("total_GERP_POS", "total_aligned_POS"))
con$variable <- as.character(con$variable)
con$variable[con$variable == "total_GERP_POS"] <- "GERP"
con$variable[con$variable == "total_aligned_POS"] <- "Aligned"


# calculate fraction 
con$fraction <- con$value/con$total_percentile_POS
con$fold_change <- con$value/con$total_percentile_POS

# group percentiles with custom function
h_tmp <- ddply(con[con$species == "Human",], "variable", gobbler)
m_tmp <- ddply(con[con$species == "Mouse",], "variable", gobbler)
dt_plot <- rbind(h_tmp, m_tmp)

# set levels order
levels(dt_plot$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                                      "51-75", "76-90", "91-95", "96-98", "99", "100")

# create group
dt_plot$group[dt_plot$variable == "GERP" & dt_plot$species == "Human"] <- "Human GERP" 
dt_plot$group[dt_plot$variable == "Aligned" & dt_plot$species == "Human"] <- "Human Alignment" 
dt_plot$group[dt_plot$variable == "GERP" & dt_plot$species == "Mouse"] <- "Mouse GERP" 
dt_plot$group[dt_plot$variable == "Aligned" & dt_plot$species == "Mouse"] <- "Mouse Alignment" 

# # calculate fold change vs the 100th percentile
# con$fold_change[con$species == "Human" & con$variable == "GERP"] <- con$fraction[con$species == "Human" & con$variable == "GERP"] / 
#                                                 con$fraction[con$species == "Human" & con$variable == "GERP" & con$percentile == 100]
# con$fold_change[con$species == "Human" & con$variable == "Aligned"] <- con$fraction[con$species == "Human" & con$variable == "Aligned"] / 
#                                                 con$fraction[con$species == "Human" & con$variable == "Aligned" & con$percentile == 100]
# con$fold_change[con$species == "Mouse" & con$variable == "GERP"] <- con$fraction[con$species == "Human" & con$variable == "GERP"] / 
#                                                 con$fraction[con$species == "Mouse" & con$variable == "GERP" & con$percentile == 100]
# con$fold_change[con$species == "Mouse" & con$variable == "Aligned"] <- con$fraction[con$species == "Mouse" & con$variable == "Aligned"] / 
#                                                 con$fraction[con$species == "Mouse" & con$variable == "Aligned" & con$percentile == 100]
# 
# # group percentiles with custom function
# h_tmp <- ddply(con[con$species == "Human",], "variable", gobbler)
# m_tmp <- ddply(con[con$species == "Mouse",], "variable", gobbler)
# plot_dt <- rbind(h_tmp, m_tmp)
# 
# # set levels order
# levels(plot_dt$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                                        "51-75", "76-90", "91-95", "96-98", "99", "100")

### PLOT by fraction

p1 <- ggplot(dt_plot, aes(x=rank, y=mean_fold_change, color=group)) + 
  # geom_point() +
  geom_line(size=1.5) +
  geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change), 
                stat = "identity",
                width=0.2,
                size=1) + 
  xlab("Constraint (percentile rank)") +
  ylab("Territory fraction") +
  # scale_y_continuous(trans = log2_trans(),
  #                    # breaks = trans_breaks("log2", function(x) 2^x)) +
  #                    breaks = c(1, 2, 5, 10, 20, 40),
  #                    limits = c(1, 60)) +
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
p1
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_conservation_by_percentile.jpg", plot = p1, height = 5, width = 6)


####

# ### TEST ALIGNMENT
# 
# M_con <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_mouse_950_50.csv")
# # get alignment fraction for each percentile
# M_con$length <- (M_con$POS_to + 1) - M_con$POS_from
# M_frac <- merge(
#   aggregate(n_aligned_central ~ Constraint_percentile_CpG, data=M_con, FUN=sum),
#   aggregate(length ~ Constraint_percentile_CpG, data=M_con, FUN=sum)
# )
# M_frac$fraction <- M_frac$n_aligned_central / M_frac$length

# dt_plot <- subset(dt_plot, dt_plot$group == "Human Alignment" | dt_plot$group == "Mouse Alignment")
# 
# p2 <- ggplot(dt_plot, aes(x=rank, y=mean_fold_change, color=species)) + 
#   # geom_point() +
#   geom_line(size=1.5) +
#   geom_errorbar(aes(ymax = upper_range_fold_change, ymin = lower_range_fold_change), 
#                 stat = "identity",
#                 width=0.2,
#                 size=1) + 
#   xlab("Constraint (percentile rank)") +
#   ylab("Human-mouse conservation fraction") +
#   # scale_y_continuous(trans = log2_trans(),
#   #                    # breaks = trans_breaks("log2", function(x) 2^x)) +
#   #                    breaks = c(1, 2, 5, 10, 20, 40),
#   #                    limits = c(1, 60)) +
#   scale_x_continuous(breaks=c(1:12),
#                      labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
#                               "51-75", "76-90", "91-95", "96-98", "99", "100")) +
#   theme_bw() +
#   theme(
#     legend.title = element_blank(),
#     legend.justification=c(1,0), 
#     legend.position=c(0.95, 0.70),  
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     text = element_text(size=14),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# p2
# ggsave("~/Dropbox/Figure_conservation_by_percentile.jpg", plot = p2, height = 5, width = 6)
# 

