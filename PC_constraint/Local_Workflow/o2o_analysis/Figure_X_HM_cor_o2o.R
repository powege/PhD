### FIGURE 2 --

### Add regression equation and r2

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)

### FUNCTION for plotting lm equation on ggplot
linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 2),
            yy = format(abs(coef(k)[2]), digits = 2),
            r2 = format(summary(k)$r.squared, digits = 3));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)   
  }
  as.character(as.expression(eq));               
}

### FUNCTION for plotting Pearson's correlation on ggplot
corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}

### IMPORT data
df <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")

### FORMAT 

# caclulate percentile for orhtologues
percentile <- ecdf(df$M_fun_Z_0.0001[!duplicated(df$M_external_gene_name)]) # remove duplicates 
df$M_fun_Z_percentile_orth <- percentile(df$M_fun_Z_0.0001)
# hist(df$M_fun_Z_percentile_orth)
percentile <- ecdf(df$H_fun_Z_0.001[!duplicated(df$H_external_gene_name)]) # remove duplicates
df$H_fun_Z_percentile_orth <- percentile(df$H_fun_Z_0.001)
# hist(df$H_fun_Z_percentile_orth)

df.plot <- df[!is.na(df$M_fun_Z_0.0001) & !is.na(df$H_fun_Z_0.001),]
df.plot$HM_funZ_class <- as.factor(df.plot$HM_funZ_class)
df.plot$HM_funZ_class <- factor(df.plot$HM_funZ_class, levels = c("NA", "Most constrained", "Least constrained"))
summary(df.plot)


### MODEL

mod <- lm(df.plot$H_fun_Z_0.001 ~ df.plot$M_fun_Z_0.0001)
# summary(mod)
# plot(mod)
cor.test(df.plot$H_fun_Z_0.001, df.plot$M_fun_Z_0.0001)

pearsons <- corr_eqn(df.plot$H_fun_Z_0.001, df.plot$M_fun_Z_0.0001, digits = 3)
lm_equation <- linear(mod)

df.plot$orthology_type <- factor(df.plot$orthology_type, levels = c("ortholog_one2one", "ortholog_one2many", "ortholog_many2many"))

Fig2 <- ggplot(df.plot, aes(x = M_fun_Z_0.0001, y = H_fun_Z_0.001)) +
  geom_point(aes(col=orthology_type, alpha=orthology_type)) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
  # annotate("text", x = -4, y = 8, label = lm_equation, colour="black", size = 5, parse=TRUE) +
  annotate("text", x = -7, y = 8, label = pearsons, colour="black", size = 5, parse=TRUE) +
  scale_color_manual(breaks = c("ortholog_one2one", "ortholog_one2many", "ortholog_many2many"), values=c('gray80', 'red', 'blue')) +
  scale_alpha_manual(values=c(0.3,0.8,0.8)) + 
  xlab("Mouse funZ") +
  ylab('Human funZ') +
  xlim(-10, 10) +
  ylim(-10, 10) +
  theme_bw() +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
Fig2

ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_X_HM_cor_orth_type.jpg", plot = Fig2, width = 6, height = 6)
cor <- cor.test(df.plot$H_fun_Z_0.001, df.plot$M_fun_Z_0.0001, method = "pearson", exact = T)
cor
cor$p.value

summary(df.plot)
HMC <- subset(df, df$HM_funZ_class == "Most constrained")
length(unique(HMC$H_external_gene_name))
length(unique(HMC$M_external_gene_name))
HMU <- subset(df, df$HM_funZ_class == "Least constrained")
length(unique(HMU$H_external_gene_name))
length(unique(HMU$M_external_gene_name))

######################

### FORMAT 

# subset 1to1 orthologs
df <- subset(df, df$orthology_type == "ortholog_one2one")

# caclulate percentile for orhtologues
percentile <- ecdf(df$M_fun_Z_0.0001[!duplicated(df$M_external_gene_name)]) # remove duplicates 
df$M_fun_Z_percentile_orth <- percentile(df$M_fun_Z_0.0001)
# hist(df$M_fun_Z_percentile_orth)
percentile <- ecdf(df$H_fun_Z_0.001[!duplicated(df$H_external_gene_name)]) # remove duplicates
df$H_fun_Z_percentile_orth <- percentile(df$H_fun_Z_0.001)
# hist(df$H_fun_Z_percentile_orth)

### identify overlapping extremes as most tolerant and intolerant
df$HM_funZ_class <- "NA"
df$HM_funZ_class[df$H_fun_Z_percentile_0.001 >= 0.9 & df$M_fun_Z_percentile_0.0001 >= 0.9] <- "Most constrained"
df$HM_funZ_class[df$H_fun_Z_percentile_0.001 <= 0.1 & df$M_fun_Z_percentile_0.0001 <= 0.1] <- "Least constrained"
# df$HM_funZ_class[df$H_fun_Z_percentile_orth >= 0.9 & df$M_fun_Z_percentile_orth >= 0.9] <- "Human and mouse constrained"
# df$HM_funZ_class[df$H_fun_Z_percentile_orth <= 0.1 & df$M_fun_Z_percentile_orth <= 0.1] <- "Human and mouse unconstrained"


# remove extreme values
# M_removed <- subset(df, df$M_fun_Z_0.0001 > 10 |  df$M_fun_Z_0.0001 < -10 )
# H_removed <- subset(df, df$H_fun_Z_0.001 > 10 |  df$H_fun_Z_0.001 < -10 )

df.plot <- df[!is.na(df$M_fun_Z_0.0001) & !is.na(df$H_fun_Z_0.001),]
df.plot$HM_funZ_class <- as.factor(df.plot$HM_funZ_class)
df.plot$HM_funZ_class <- factor(df.plot$HM_funZ_class, levels = c("NA", "Most constrained", "Least constrained"))
summary(df.plot)


### MODEL

mod <- lm(df.plot$H_fun_Z_0.001 ~ df.plot$M_fun_Z_0.0001)
# summary(mod)
# plot(mod)
cor.test(df.plot$H_fun_Z_0.001, df.plot$M_fun_Z_0.0001)

pearsons <- corr_eqn(df.plot$H_fun_Z_0.001, df.plot$M_fun_Z_0.0001, digits = 3)
lm_equation <- linear(mod)

Fig2 <- ggplot(df.plot, aes(x = M_fun_Z_0.0001, y = H_fun_Z_0.001)) +
  geom_point(aes(col=HM_funZ_class)) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
  # annotate("text", x = -4, y = 8, label = lm_equation, colour="black", size = 5, parse=TRUE) +
  annotate("text", x = -7, y = 8, label = pearsons, colour="black", size = 5, parse=TRUE) +
  scale_color_manual(breaks = c("Most constrained", "Least constrained"), values=c('gray80', 'red', 'blue')) +
  xlab("Mouse funZ") +
  ylab('Human funZ') +
  xlim(-10, 10) +
  ylim(-10, 10) +
  theme_bw() +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
Fig2

# ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_2_o2o.jpg", plot = Fig2, width = 6, height = 6)
cor <- cor.test(df.plot$H_fun_Z_0.001, df.plot$M_fun_Z_0.0001, method = "pearson", exact = T)
cor
cor$p.value

summary(df.plot)
HMC <- subset(df, df$HM_funZ_class == "Most constrained")
length(unique(HMC$H_external_gene_name))
length(unique(HMC$M_external_gene_name))
HMU <- subset(df, df$HM_funZ_class == "Least constrained")
length(unique(HMU$H_external_gene_name))
length(unique(HMU$M_external_gene_name))



