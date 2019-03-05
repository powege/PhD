### SCRIPT that analyses ClinVar variant counts and gene constraint

rm(list = ls())
graphics.off()

library(data.table)
library("cowplot")

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

### FUNCTION that plots mean and median ClinVar pathogenic variants agaisnt funZ percentile
M_plotter <- function(m_df_md, m_df_mn, lm_equation){
  
  m_plot <- ggplot() +
    geom_bar(data=m_df_md, aes(y=x,x=Category), color = "red", fill="red", stat="identity", width = 0.75) + 
    geom_point(data=m_df_mn, aes(y=x, x=Category), alpha = 1/5) +
    geom_smooth(data=m_df_mn, aes(y=x, x=Category), color = "blue", method='lm') +
    # scale_color_gradient(low="blue", high="red") +
    annotate("text", x = 39, y = 3.5, label = lm_equation, colour="black", size = 5, parse=TRUE) +
    xlab("Mouse funZ percentile bin") +
    ylab('Pathogenic variants per kb\nin human orthologue') +ylim(0, 4) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 14),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          plot.margin=unit(c(1.5,1,1,1),"cm"))
  
  return(m_plot)
}

### FUNCTION that plots mean and median ClinVar pathogenic variants agaisnt funZ percentile
H_plotter <- function(m_df_md, m_df_mn, lm_equation){
  
  m_plot <- ggplot() +
    geom_bar(data=m_df_md, aes(y=x,x=Category), color = "red", fill="red", stat="identity", width = 0.75) + 
    geom_point(data=m_df_mn, aes(y=x, x=Category), alpha = 1/5) +
    geom_smooth(data=m_df_mn, aes(y=x, x=Category), color = "blue", method='lm') +
    # scale_color_gradient(low="blue", high="red") +
    annotate("text", x = 39, y = 3.5, label = lm_equation, colour="black", size = 5, parse=TRUE) +
    xlab("Human funZ percentile bin") +
    ylab('Pathogenic variants per kb') +ylim(0, 4) +
    theme_bw() +
    theme(legend.position="none",
          text = element_text(size = 14),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          plot.margin=unit(c(1.5,1,1,1),"cm"))
  
  return(m_plot)
}


# IMPORT
CV <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/ClinVar_VEP_output.txt", fill = T)
con <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
len <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_funZ_0.001.csv")

### FORMAT

con <- con[con$orthology_type == "ortholog_one2one"]

colnames(CV) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "Gene",
                  "Feature", "Feature_type", "Consequence", "IMPACT", "SYMBOL",
                  "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", 'CCDS')
CV <- CV[CV$Consequence %like% "synonymous_variant" |
           CV$Consequence %like% "stop_retained_variant" |
           CV$Consequence %like% "start_retained_variant" |
           CV$Consequence %like% "missense_variant" |
           CV$Consequence %like% "stop_gained" |
           CV$Consequence %like% "start_lost" |
           CV$Consequence %like% "stop_lost" |
           CV$Consequence %like% "splice_donor_variant" |
           CV$Consequence %like% "splice_acceptor_variant",]
CV <- as.data.frame(table(CV$SYMBOL))
colnames(CV) <- c("H_external_gene_name", "n_ClinVar")

# number of variants
sum(CV$n_ClinVar)

# subset cds length
len <- len[,c("external_gene_name", "cds_length")]
colnames(len) <- paste0("H_", colnames(len))

# merge
df <- merge(con, len)
df <- merge(df, CV, all = T)

# subset one 2 one ortholgs
# df <- df[df$orthology_type == "ortholog_one2one",]

# subset mouse genes
m_df <- df[!duplicated(df$M_external_gene_name),]

# calculate one to one percentile
percentile <- ecdf(m_df$M_fun_Z_0.0001)
m_df$M_fun_Z_percentile <- percentile(m_df$M_fun_Z_0.0001)
m_df$M_fun_Z_percentile <- ceiling((m_df$M_fun_Z_percentile)*100)

# calculate n_ClinVar per kb (n_ClinVar/cds_length)*1000
m_df$rel_ClinVar <- (m_df$n_ClinVar/m_df$H_cds_length)*1000
# length(unique(m_df$M_external_gene_name))

# calculate rel_ClinVar per gene
m_df <- m_df[,c("M_fun_Z_percentile", "rel_ClinVar")]
m_df$rel_ClinVar[is.na(m_df$rel_ClinVar)] <- 0

# m_df <- aggregate(m_df$rel_ClinVar, by=list(Category=m_df$M_fun_Z_percentile), FUN=sum)
m_df_mn <- aggregate(m_df$rel_ClinVar, by=list(Category=m_df$M_fun_Z_percentile), FUN=mean)
m_df_md <- aggregate(m_df$rel_ClinVar, by=list(Category=m_df$M_fun_Z_percentile), FUN=median)

# fit model
m_mn_mod <- lm(x ~ Category, data = m_df_mn)
summary(m_mn_mod)
m_lm_equation <- linear(m_mn_mod)

# subset human genes
h_df <- df[!duplicated(df$H_external_gene_name),]

# calculate one to one percentile
percentile <- ecdf(h_df$M_fun_Z_0.0001)
h_df$M_fun_Z_percentile <- percentile(h_df$M_fun_Z_0.0001)
h_df$M_fun_Z_percentile <- ceiling((h_df$M_fun_Z_percentile)*100)

# calculate n_ClinVar per kb (n_ClinVar/cds_length)*1000
h_df$rel_ClinVar <- (h_df$n_ClinVar/h_df$H_cds_length)*1000
# length(unique(h_df$H_external_gene_name))

# calculate rel_ClinVar per gene
h_df <- h_df[,c("M_fun_Z_percentile", "rel_ClinVar")]
h_df$rel_ClinVar[is.na(h_df$rel_ClinVar)] <- 0

# h_df <- aggregate(h_df$rel_ClinVar, by=list(Category=h_df$M_fun_Z_percentile), FUN=sum)
h_df_md <- aggregate(h_df$rel_ClinVar, by=list(Category=h_df$M_fun_Z_percentile), FUN=median)
h_df_mn <- aggregate(h_df$rel_ClinVar, by=list(Category=h_df$M_fun_Z_percentile), FUN=mean)

# fit model
h_mn_mod <- lm(x ~ Category, data = h_df_mn)
summary(h_mn_mod)
h_lm_equation <- linear(h_mn_mod)
# plot(mod)

m_plot <- M_plotter(m_df_md, m_df_mn, m_lm_equation)
m_plot
h_plot <- H_plotter(h_df_md, h_df_mn, h_lm_equation)
h_plot

Fig <- plot_grid(m_plot, h_plot, ncol = 2, nrow = 1, labels = c("A", "B"), 
                 label_size = 20, hjust = -0.1, vjust = 2)
Fig
save_plot("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_ClinVar_o2o.jpg", Fig, ncol = 1, nrow = 1, base_height = 5, base_width = 10)

