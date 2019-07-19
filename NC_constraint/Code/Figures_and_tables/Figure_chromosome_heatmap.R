rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)
library("cowplot")
library(gridExtra)


### HUMAN
h_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_human_950_50_MAF001.csv")
# dummy <- subset(dt, dt$CHR %in% c(21, 22))

# format
h_dt$POS_from <- h_dt$POS_from / 1000000
h_dt$POS_to <- h_dt$POS_to / 1000000

p_human <- ggplot(h_dt) + 
  geom_rect(aes(xmin = CHR - 0.3, xmax = CHR + 0.3, ymin = POS_from, ymax = POS_to, fill = Constraint_percentile_CpG)) + 
  scale_fill_distiller(palette = "RdYlBu",
                       direction = -1,
                       breaks = c(1, 25, 50, 75, 100),
                       limits = c(1, 100)) +
  # scale_y_reverse() + 
  # scale_x_reverse(breaks = c(1:22),
  #                 limits = c(0.5, 22.5)) +
  xlab("Chromosome") +
  ylab("Mb") +
  ggtitle("Human") +
  scale_x_continuous(breaks = c(1:22),
                     limits = c(0.5, 22.5) #,
                     # position = "top"
                     ) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250),
                     limits = c(0, 250)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.justification=c(1,0),
    legend.position=c(0.95, 0.65),
    # legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  ) 
# p_human
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_constraint_heat_map_human.jpg", plot = p_human, height = 5, width = 6)


### MOUSE
m_dt <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_window_mouse_950_50.csv")
# dummy <- subset(dt, dt$CHR %in% c(21, 22))

# format
m_dt$POS_from <- m_dt$POS_from / 1000000
m_dt$POS_to <- m_dt$POS_to / 1000000

p_mouse <- ggplot(m_dt) + 
  geom_rect(aes(xmin = CHR - 0.3, xmax = CHR + 0.3, ymin = POS_from, ymax = POS_to, fill = Constraint_percentile_CpG)) + 
  # scale_y_reverse() + 
  scale_fill_distiller(palette = "RdYlBu",
                       direction = -1,
                       breaks = c(1, 25, 50, 75, 100),
                       limits = c(1, 100)) +
  xlab("Chromosome") +
  ylab("Mb") +
  ggtitle("Mouse") +
  scale_x_continuous(breaks = c(1:19),
                     limits = c(0.5, 19.5) #,
                     # position = "top"
                     ) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250),
                     limits = c(0, 250)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.justification=c(1,0), 
    legend.position=c(0.95, 0.65), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  ) 
# p_mouse
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_constraint_heat_map_mouse.jpg", plot = p_mouse, height = 5, width = 6)


pout <- grid.arrange(p_mouse, p_human, nrow = 1, widths = c(6.5, 6.5))
ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_HM_SNV_constraint_heat_map.jpg", plot = pout, height = 6, width = 10)









