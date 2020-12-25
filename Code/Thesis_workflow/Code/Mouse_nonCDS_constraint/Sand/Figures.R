rm(list=ls())
graphics.off()

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)


tss <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Sand/m.mus_GRC38_Ensembl_v101_regulatory_features_multicell_closest_TSS.csv.gz")
con <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Sand/m.mus_GRC38_Ensembl_v101_regulatory_features_multicell_constraint.csv.gz")
impc <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/IMPC/Formatted/IMPC_r12_viability_fertility_pleiotropy.csv")
ens <- fread("~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/BioMart/Ensembl_v94_m.mus_canPC_pos.csv")

# colnames(tss)
# colnames(con)
# table(tss$category)
# table(con$category)

ens <- unique(ens[,c("external_gene_name", "ensembl_transcript_id")])

colnames(tss) <- c("chromosome",
                   "start",
                   "end",
                   "category",
                   "category_id",
                   "ensembl_transcript_id",
                   "transcription_start_site",
                   "TSS_dist",
                   "domain",
                   "f_category_overlap",
                   "f_domain_overlap")

con <- con[,c("chromosome",
              "start",
              "end",
              "category",
              "category_id",
              "transcript_id",
              "percentage_Nmask",            
              "percentage_sm",
              "n_snv",
              "n_Nmask",
              "n_DPmask",
              "n_soft_mask",
              "n_CG",
              "n_meth_CG",
              "length",
              "percentage_DPmask",           
              "exp_snv",
              "exp_snv_lwr",
              "exp_snv_upr",                 
              "OE_snv",
              "OE_snv_lwr",
              "OE_snv_upr")]
colnames(con)[names(con) == "transcript_id"] <- "ensembl_transcript_id"
con$ensembl_transcript_id[con$ensembl_transcript_id == ""] <- NA

con$percent_CG <- ( con$n_CG / con$length ) * 100
con$percent_CG_meth <- ( con$n_meth_CG / con$n_CG ) * 100

dt <- full_join(tss, con, by = c("chromosome",
                                 "start",
                                 "end",
                                 "category",
                                 "category_id"))

dt$ensembl_transcript_id <- NA
dt$ensembl_transcript_id[dt$category != "Exon - UTR"] <- dt$ensembl_transcript_id.x[dt$category != "Exon - UTR"]
dt$ensembl_transcript_id[dt$category == "Exon - UTR"] <- dt$ensembl_transcript_id.y[dt$category == "Exon - UTR"]
dt$ensembl_transcript_id.x <- NULL
dt$ensembl_transcript_id.y <- NULL
dt <- as.data.table(dt)
dt <- dt[ens, on = "ensembl_transcript_id"]

dt <- dt[impc, on = c("external_gene_name")]
dt <- dt[!is.na(category)]
dt <- dt[!is.na(ensembl_transcript_id)]

out_list <- list()
for (i in 1:length(unique(dt$category))){
  
  sub <- dt[category == unique(dt$category)[i]]
  
  percentile <- ecdf(sub$OE_snv_lwr[!duplicated(sub$external_gene_name)])
  sub$OE_lwr_percentile <- ceiling(percentile(sub$OE_snv_lwr)*100)
  
  percentile <- ecdf(sub$OE_snv[!duplicated(sub$external_gene_name)])
  sub$OE_percentile <- ceiling(percentile(sub$OE_snv)*100)
  
  percentile <- ecdf(sub$OE_snv_upr[!duplicated(sub$external_gene_name)])
  sub$OE_upr_percentile <- ceiling(percentile(sub$OE_snv_upr)*100)
  
  # percentile <- ecdf(sub$RVIS_specific[!duplicated(sub$external_gene_name)])
  # sub$RVIS_percentile <- ceiling(percentile(sub$RVIS_specific)*100)
  
  out_list[[i]] <- sub
}
dt <- do.call("rbind", out_list)

# set category of interest
dt_plot <- dt[category == "Enhancer - proximal"]
dt_plot <- dt[category == "Enhancer - distal"]
dt_plot <- dt[category == "Exon - UTR"]
dt_plot <- dt[category == "Promoter"]


### viability status

hist(dt_plot$n_CG)
my_comparisons <- list( c("L", "VN") )
ggplot(dt_plot, aes(x=viability_status, y=n_CG, fill=viability_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

hist(dt_plot$percent_CG)
ggplot(dt_plot, aes(x=viability_status, y=percent_CG, fill=viability_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

hist(dt_plot$percent_CG_meth)
ggplot(dt_plot, aes(x=viability_status, y=percent_CG_meth, fill=viability_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

hist(dt_plot$OE_percentile)
ggplot(dt_plot, aes(x=viability_status, y=OE_percentile, fill=viability_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")


### fertility status

length(unique(dt$ensembl_transcript_id[dt$fertility_status == "FFMF"]))
length(unique(dt$ensembl_transcript_id[dt$fertility_status == "FI"]))
length(unique(dt$ensembl_transcript_id[dt$fertility_status == "FIMI"]))
length(unique(dt$ensembl_transcript_id[dt$fertility_status == "MI"]))

hist(dt_plot$n_CG)
my_comparisons <- list( c("FFMF", "FIMI"), c("FFMF", "FI"), c("FFMF", "MI") )
ggplot(subset(dt_plot, fertility_status != ""), aes(x=fertility_status, y=n_CG, fill=fertility_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

hist(dt_plot$percent_CG)
ggplot(subset(dt_plot, fertility_status != ""), aes(x=fertility_status, y=percent_CG, fill=fertility_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

hist(dt_plot$percent_CG_meth)
ggplot(subset(dt_plot, fertility_status != ""), aes(x=fertility_status, y=percent_CG_meth, fill=fertility_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

hist(dt_plot$OE_percentile)
ggplot(subset(dt_plot, fertility_status != ""), aes(x=fertility_status, y=OE_percentile, fill=fertility_status)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")


## n phenotypes

dt_plot_mn <- aggregate(dt_plot$hit_rate_procedure, by=list(Category=dt_plot$OE_percentile), FUN=mean, na.rm = TRUE)
cor.test(dt_plot_mn$Category, dt_plot_mn$x)
ggplot(dt_plot_mn, aes(x=Category, y=x)) +
  geom_point() + 
  geom_smooth(method = "lm")

dt_plot_mn <- aggregate(dt_plot$hit_rate_parameter, by=list(Category=dt_plot$OE_percentile), FUN=mean, na.rm = TRUE)
cor.test(dt_plot_mn$Category, dt_plot_mn$x)
ggplot(dt_plot_mn, aes(x=Category, y=x)) +
  geom_point() + 
  geom_smooth(method = "lm")

hist(dt_plot$percent_CG_meth)
dt_plot$percent_CG_meth_bin <- NA
dt_plot$percent_CG_meth_bin[dt_plot$percent_CG_meth < 60 ] <- "Low"
dt_plot$percent_CG_meth_bin[dt_plot$percent_CG_meth >= 60 ] <- "High"

hist(dt_plot$hit_rate_parameter)
cor.test(dt_plot$hit_rate_parameter, dt_plot$percent_CG_meth, method = "spearman")
ggplot(subset(dt_plot, !is.na(percent_CG_meth_bin)), aes(x=percent_CG_meth_bin, y=hit_rate_parameter)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

hist(dt_plot$hit_rate_procedure)
cor.test(dt_plot$hit_rate_procedure, dt_plot$percent_CG_meth, method = "spearman")
ggplot(subset(dt_plot, !is.na(percent_CG_meth_bin)), aes(x=percent_CG_meth_bin, y=hit_rate_procedure, fill=percent_CG_meth_bin)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

hist(dt_plot$hit_rate_top_level_mp)
cor.test(dt_plot$hit_rate_top_level_mp, dt_plot$percent_CG_meth, method = "spearman")
ggplot(subset(dt_plot, !is.na(percent_CG_meth_bin)), aes(x=percent_CG_meth_bin, y=hit_rate_top_level_mp)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")










