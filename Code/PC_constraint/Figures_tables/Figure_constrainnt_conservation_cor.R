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
spear_corr_eqn <- function(x, y, digits = 2) {
  corr_coef <- round(cor(x, y, method = "spearman"), digits = digits)
  paste("italic(rho) == ", corr_coef)
}

### IMPORT data
df <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/HM_constriant_orth_v2.csv")
orths <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv")


### FORMAT 

orths <- orths[,c("H_ensembl_transcript_id", "Maa_match_Haa", "Haa_match_Maa", "M_ensembl_transcript_id")]
orths <- transform(orths, Conservation = pmin(Haa_match_Maa, Maa_match_Haa))
# orths$Conservation <- (orths$Maa_match_Haa + orths$Haa_match_Maa)/2

# merge
df <- orths[df, on = c("H_ensembl_transcript_id", "M_ensembl_transcript_id")]

# subset 1to1 orthologs
# df <- subset(df, df$orthology_type == "ortholog_one2one")

# caclulate percentile for orhtologues
percentile <- ecdf(df$M_fun_Z[!duplicated(df$M_external_gene_name)]) # remove duplicates 
df$M_fun_Z_percentile <- percentile(df$M_fun_Z)
# hist(df$M_fun_Z_percentile_orth)
percentile <- ecdf(df$H_fun_Z[!duplicated(df$H_external_gene_name)]) # remove duplicates
df$H_fun_Z_percentile <- percentile(df$H_fun_Z)
# hist(df$H_fun_Z_percentile_orth)

### identify overlapping extremes as most tolerant and intolerant
df$HM_funZ_class <- "NA"
df$HM_funZ_class[df$H_fun_Z_percentile >= 0.9 & df$M_fun_Z_percentile >= 0.9] <- "Most constrained"
df$HM_funZ_class[df$H_fun_Z_percentile <= 0.1 & df$M_fun_Z_percentile <= 0.1] <- "Least constrained"
# df$HM_funZ_class[df$H_fun_Z_percentile_orth >= 0.9 & df$M_fun_Z_percentile_orth >= 0.9] <- "Human and mouse constrained"
# df$HM_funZ_class[df$H_fun_Z_percentile_orth <= 0.1 & df$M_fun_Z_percentile_orth <= 0.1] <- "Human and mouse unconstrained"


# remove extreme values
# M_removed <- subset(df, df$M_fun_Z_0.0001 > 10 |  df$M_fun_Z_0.0001 < -10 )
# H_removed <- subset(df, df$H_fun_Z_0.001 > 10 |  df$H_fun_Z_0.001 < -10 )

df.plot <- df[!is.na(df$M_fun_Z) & !is.na(df$H_fun_Z),]
df.plot$HM_funZ_class <- as.factor(df.plot$HM_funZ_class)
df.plot$HM_funZ_class <- factor(df.plot$HM_funZ_class, levels = c("NA", "Most constrained", "Least constrained"))
summary(df.plot)


### MODEL

cor.test(df.plot$Conservation, df.plot$M_fun_Z, method = "spearman")
spearman <- spear_corr_eqn(df.plot$Conservation, df.plot$M_fun_Z, digits = 3)

plot(df.plot$Conservation, df.plot$M_OE_ratio)
# cor.test(df.plot$Conservation, df.plot$M_OE_ratio, method = "spearman")


M_Fig <- ggplot(df.plot, aes(x = M_fun_Z, y = Conservation)) +
  geom_point(aes(col=HM_funZ_class)) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
  # annotate("text", x = -4, y = 8, label = lm_equation, colour="black", size = 5, parse=TRUE) +
  annotate("text", x = -7, y = 8, label = spearman, colour="black", size = 5, parse=TRUE) +
  scale_color_manual(breaks = c("Most constrained", "Least constrained"), values=c('gray80', 'red', 'blue')) +
  xlab("Mouse funZ") +
  ylab('Conservation (%)') +
  xlim(-10, 10) +
  ylim(0, 100) +
  theme_bw() +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
M_Fig

cor.test(df.plot$Conservation, df.plot$H_fun_Z, method = "spearman")
spearman <- spear_corr_eqn(df.plot$Conservation, df.plot$H_fun_Z, digits = 3)
# cor.test(df.plot$Conservation, df.plot$H_OE_ratio, method = "spearman")

H_Fig <- ggplot(df.plot, aes(x = H_fun_Z, y = Conservation)) +
  geom_point(aes(col=HM_funZ_class)) +
  geom_smooth(method='lm', formula=y~x, se = T, colour = "black", size = 0.6, fullrange = T) +
  # annotate("text", x = -4, y = 8, label = lm_equation, colour="black", size = 5, parse=TRUE) +
  annotate("text", x = -7, y = 8, label = spearman, colour="black", size = 5, parse=TRUE) +
  scale_color_manual(breaks = c("Most constrained", "Least constrained"), values=c('gray80', 'red', 'blue')) +
  xlab("Human funZ") +
  ylab('Conservation (%)') +
  xlim(-10, 10) +
  ylim(0, 100) +
  theme_bw() +
  theme(legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        text = element_text(size = 14),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"))
H_Fig

## Combine
out <- plot_grid(M_Fig, H_Fig, labels = c("Mouse", "Human"),
                 ncol = 2, nrow = 1)
out
save_plot("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/PLoS/Figure_funZ_conservation_cor.jpg", 
          out, ncol = 2, nrow = 1, base_height = 5)



