rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)


### IMPORT

h_dt <- list()
for (i in 1:22){
  h_dt[[i]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_CpG_proportion_human_chr", i, ".csv"))
}
h_dt <- do.call("rbind", h_dt)

m_dt <- list()
for (i in 1:19){
  m_dt[[i]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/Figure_annotation_CpG_proportion_mouse_chr", i, ".csv"))
}
m_dt <- do.call("rbind", m_dt)

### Format

# add species
h_dt$Species <- "Human"
m_dt$Species <- "Mouse"

# rbind
dt <- rbind(h_dt, m_dt)
rm(h_dt, m_dt)

# QC
plot_dt <- subset(dt, dt$N_proportion == 0 & dt$length >= 4)

# calculate 1st 2nd and 3rd quartiles for each annotation
tmp1 <- aggregate(plot_dt[, 3], list(plot_dt$Species, plot_dt$category), median)
colnames(tmp1) <- c("Species", "Category", "Q2")
tmp2 <- aggregate(plot_dt[, 3], list(plot_dt$Species, plot_dt$category), quantile, probs = 0.25)
colnames(tmp2) <- c("Species", "Category", "Q1")
tmp3 <- aggregate(plot_dt[, 3], list(plot_dt$Species, plot_dt$category), quantile, probs = 0.75)
colnames(tmp3) <- c("Species", "Category", "Q3")
plot_tmp <- merge(tmp2, tmp1)
plot_tmp <- merge(plot_tmp, tmp3)

# set factor order by median CpG fraction
order <- aggregate(plot_tmp$Q2~plot_tmp$Category, FUN=median)
order <- order[order(order$`plot_tmp$Q2`),]
plot_tmp$Category <- factor(plot_tmp$Category, levels = as.character(order$`plot_tmp$Category`))
# plot_tmp$Species <- factor(plot_tmp$Species, levels = c("Human", "Mouse"))


# set rectangle coordinates 
rects <- data.frame(xmin = head(seq <- seq(0.5, 10 + 0.5, 1), -1), 
                    xmax = tail(seq, -1), 
                    Category = levels(plot_tmp$Category),
                    rect_type = c("a", "c"))
rects <- rbind(rects)
plot_tmp <- merge(plot_tmp, rects)


### PLOT

p1 <- ggplot(plot_tmp, aes(x=Category, y=Q2, color=Species)) + 
  geom_errorbar(aes(ymax = Q3, ymin = Q1), 
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.5,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
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
                position = position_dodge(width=0.9), 
                stat = "identity",
                width=0.4,
                size=1.2) + 
  geom_point(position = position_dodge(0.9), size = 2.5) + 
  scale_fill_manual(values = c("grey", "white")) +
  xlab("Genomic annotation") +
  ylab("CpG fraction") +
  scale_y_continuous(breaks = c(0, 0.04, 0.08, 0.12, 0.16),
                     limits = c(0, 0.165)) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
p1

ggsave("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_CpG_composition.jpg", plot = p1, height = 5, width = 7)


#####

# p <- ggplot(plot_dt, aes(x=category, y=CpG_proportion, fill=Species)) +
#   geom_boxplot() + 
#   coord_flip()
# p

# m_prom <- subset(dt, dt$Species == "Mouse" & dt$category == "Promoter")
# h_prom <- subset(dt, dt$Species == "Human" & dt$category == "Promoter")
# mean(m_prom$length)
# mean(h_prom$length)

