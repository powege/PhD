rm(list=ls())
graphics.off()

library(ggplot2)

df <- read.csv("~/Dropbox/PhD/Data/NC_constraint/tmp_CDS_by_prcentile.csv")

df$Outlier <- 0
df$Outlier[df$Percentile == 1] <- 1
df$Outlier <- as.factor(df$Outlier)

ggplot(df, aes(x=Percentile, y=CDS_fraction, color=Outlier)) + 
  geom_point() +
  scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint percentile") +
  ylab("Coding sequence proportion") +
  theme(legend.position = "none")
