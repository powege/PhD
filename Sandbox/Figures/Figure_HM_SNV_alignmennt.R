rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)


### IMPORT
cv_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_z_scores_v2.csv")
gwas_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/GWAS_z_scores_v2.csv")


### FORMAT
p1_cv <- cv_z[,c("H_annotation", "F_align_z", "F_conserved_z", "F_align_ann_z")]
p1_cv <- melt(p1_cv, id.vars = "H_annotation")
p1_cv$variable <- as.character(p1_cv$variable)
p1_cv$variable[p1_cv$variable == "F_align_z"] <- "Orthologous\nbase"
p1_cv$variable[p1_cv$variable == "F_conserved_z"] <- "Conserved\nbase"
p1_cv$variable[p1_cv$variable == "F_align_ann_z"] <- "Conserved\nannotation"

p1_gwas <- gwas_z[,c("H_annotation", "F_align_z", "F_conserved_z", "F_align_ann_z")]
p1_gwas <- melt(p1_gwas, id.vars = "H_annotation")
p1_gwas$variable <- as.character(p1_gwas$variable)
p1_gwas$variable[p1_gwas$variable == "F_align_z"] <- "Orthologous\nbase"
p1_gwas$variable[p1_gwas$variable == "F_conserved_z"] <- "Conserved\nbase"
p1_gwas$variable[p1_gwas$variable == "F_align_ann_z"] <- "Conserved\nannotation"


# set factor order by total alignment
order_H <- c("Unannotated", "TF binding", "Intron", "Open chromatin", "Promoter flanking", "Enhancer",
             "Exon - non-coding", "Promoter", "Exon - UTR")
order_M <- rev(order_H)
p1_cv$H_annotation  <- factor(p1_cv$H_annotation, levels = as.character(order_H))
p1_gwas$H_annotation  <- factor(p1_gwas$H_annotation, levels = as.character(order_H))

# remove unannotated and intron
p1_cv <- subset(p1_cv, p1_cv$H_annotation != "Intron" & p1_cv$H_annotation != "Unannotated")
p1_gwas <- subset(p1_gwas, p1_gwas$H_annotation != "Intron" & p1_gwas$H_annotation != "Unannotated")


### PLOT

p_cv <- ggplot(p1_cv, aes(x=H_annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  xlab("Human genomic annotation") +
  ylab("z-score") +
  ggtitle("ClinVar") +
  coord_flip() +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20, 30, 40),
                     limits = c(-20, 41)) +
  geom_hline(yintercept=3, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=-3, linetype="dashed", color = "black", size=1) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    legend.position="bottom",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.65),
    # legend.box.background = element_rect(colour = "black"),
    # legend.spacing.x = unit(1.0, 'cm'),
    # legend.text = element_text(margin = margin(t = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) 
# + guides(fill=guide_legend(
#     keywidth=0.4,
#     keyheight=0.4,
#     default.unit="inch")
#   )
p_cv

p_gwas <- ggplot(p1_gwas, aes(x=H_annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  xlab("Human genomic annotation") +
  ylab("z-score") +
  ggtitle("GWAS") +
  coord_flip() +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20, 30, 40),
                     limits = c(-20, 41)) +
  geom_hline(yintercept=3, linetype="dashed", color = "black", size=1) +
  geom_hline(yintercept=-3, linetype="dashed", color = "black", size=1) +
  # geom_text(aes(label=value), vjust=1.6, color="black",
  #           position = position_dodge(0.9), size=3.5) +
  # annotate("text", x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2, 4.8, 5.2, 5.8, 6.2, 6.8, 7.2), 
  #          y = 30, label = c("n = 1000", "n = 100", "n = 100", "n = 100", "n = 100", "n = 100", "n = 100",
  #                            "n = 100", "n = 100", "n = 100", "n = 100", "n = 100", "n = 100", "n = 1000")) +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    legend.position="bottom",
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.65),
    # legend.box.background = element_rect(colour = "black"),
    # legend.spacing.x = unit(1.0, 'cm'),
    # legend.text = element_text(margin = margin(t = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) 
# + guides(fill=guide_legend(
#     keywidth=0.4,
#     keyheight=0.4,
#     default.unit="inch")
#   )
p_gwas

pout <- grid.arrange(p_cv, p_gwas, nrow = 1, widths = c(6, 6))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_SNV_Z_tmp.jpg", plot = pout, height = 6, width = 12)




#####
# 
# # count the number of SNVs for each annotaion:
# ann <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
# out_list <- list()
# for (i in 1:length(ann)){
#   sub <- subset(cv, cv$H_ANN %like% ann[i])
#   total_H <- nrow(sub)
#   total_M_align <- nrow(sub[!is.na(sub$M_POS),])
#   total_M_conserv <- nrow(sub[sub$H_REF == sub$M_REF,])
#   total_M_conserv_K <- nrow(sub[sub$H_Kmer == sub$M_Kmer,])
#   out_list[[i]] <-  c(  ann[i],
#                         total_H,
#                         total_M_align,
#                         total_M_conserv,
#                         total_M_conserv_K )
#   names(out_list[[i]]) <- c(  "annotation",
#                               "total_H",
#                               "total_M_align",
#                               "total_M_conserv",
#                               "total_M_conserv_K" )
# }
# 
# plot_dt1 <- as.data.frame(do.call("rbind", out_list))
# # number code annotations (plyr)
# plot_dt1$annotation <- mapvalues(plot_dt1$annotation, from=c("A","B","C","D","E","F","G","H","I","J"), 
#                               to=c("Exon - CDS",
#                                    "Exon - UTR",
#                                    "Exon - non-coding",
#                                    "Promoter",
#                                    "Enhancer",
#                                    "Open chromatin",
#                                    "TF binding",
#                                    "Promoter flanking",
#                                    "Intron",
#                                    "Unannotated"))
# cols.num <- c("total_H","total_M_align","total_M_conserv","total_M_conserv_K" )
# plot_dt1[cols.num] <- sapply(plot_dt1[cols.num],as.character)
# plot_dt1[cols.num] <- sapply(plot_dt1[cols.num],as.numeric)
# plot_dt1$frac_M_align <- plot_dt1$total_M_align/plot_dt1$total_H
# plot_dt1$frac_M_conserv <- plot_dt1$total_M_conserv/plot_dt1$total_H
# plot_dt1$frac_M_conserv_K <- plot_dt1$total_M_conserv_K/plot_dt1$total_H
# 
# # wide to long
# plot_dt1 <- melt(plot_dt1, id.vars = c("annotation","total_H","total_M_align", "total_M_conserv", "total_M_conserv_K"))
# 
# # mouse annotation composition:
# annotation <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
# H_annotation <- list()
# M_annotation <- list()
# H_ann_total <- list()
# M_ann_total <- list()
# cv_align <- cv[!is.na(cv$M_POS),]
# for (i in 1:length(annotation)){
#   sub <- subset(cv_align, cv_align$H_ANN %like% annotation[i])
#   H_ann_total[[i]] <- rep(nrow(sub), 10)
#   tmp_vec <- rep(NA, 10)
#   for (j in 1:length(annotation)){
#     tmp_vec[j] <- nrow(sub[sub$M_ANN %like% annotation[j],])
#     # print(c(i,j))
#   }
#   M_ann_total[[i]] <- tmp_vec
#   H_annotation[[i]] <- rep(annotation[i], 10)
#   M_annotation[[i]] <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
# }  
# 
# plot_dt2 <- data.table(
#                      H_annotation = unlist(H_annotation),
#                      H_ann_total = unlist(H_ann_total),
#                      M_annotation = unlist(M_annotation),
#                      M_ann_total = unlist(M_ann_total)
# )
# plot_dt2$M_ann_frac <- plot_dt2$M_ann_total/plot_dt2$H_ann_total
# 
# plot_dt2$H_annotation <- mapvalues(plot_dt2$H_annotation, from=c("A","B","C","D","E","F","G","H","I","J"),
#                                  to=c("Exon - CDS",
#                                       "Exon - UTR",
#                                       "Exon - non-coding",
#                                       "Promoter",
#                                       "Enhancer",
#                                       "Open chromatin",
#                                       "TF binding",
#                                       "Promoter flanking",
#                                       "Intron",
#                                       "Unannotated"))
# plot_dt2$M_annotation <- mapvalues(plot_dt2$M_annotation, from=c("A","B","C","D","E","F","G","H","I","J"),
#                                  to=c("Exon - CDS",
#                                       "Exon - UTR",
#                                       "Exon - non-coding",
#                                       "Promoter",
#                                       "Enhancer",
#                                       "Open chromatin",
#                                       "TF binding",
#                                       "Promoter flanking",
#                                       "Intron",
#                                       "Unannotated"))
# 
# # set factor order by total alignment
# order_H <- c("Unannotated", "TF binding", "Intron", "Open chromatin", "Promoter flanking", "Enhancer",
#              "Exon - non-coding", "Promoter", "Exon - UTR", "Exon - CDS")
# order_M <- rev(order_H)
# plot_dt2$H_annotation  <- factor(plot_dt2$H_annotation, levels = as.character(order_H))
# plot_dt2$M_annotation  <- factor(plot_dt2$M_annotation, levels = as.character(order_M))
# 
# p2 <- ggplot(data = plot_dt1, aes(x=annotation, y=value, fill=variable)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   coord_flip()
# p2
# 
# # generate null distributions for frac aligned, frac conserved, frac kmer conserved.
# # plot z scores by annotation for each variable.
# 
# 
# p1 <- ggplot() +
#   geom_bar(data=plot_dt2, aes(x=H_annotation, y=M_ann_frac, fill=M_annotation), colour = "black", stat="identity") +
#   xlab("Human genomic annotation") +
#   ylab("Mouse alignment composition (fraction)") +
#   ggtitle("B") +
#   labs(fill = "Mouse genomic\nannotation") +
#   coord_flip() +
#   scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2.0, 2.5),
#                      limits = c(0, 2.5)) +
#   scale_fill_brewer(palette="Set3") +
#   # scale_fill_brewer(palette="Paired") +
#   theme_bw() +
#   theme(
#     # legend.position = "none",
#     # legend.title = element_blank(),
#     legend.title = element_text(hjust = 0.5),
#     # legend.key.size = unit(2, 'lines'),
#     # legend.justification=c(1,0),
#     # legend.position=c(0.95, 0.05),
#     # legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(size = 26, face = "bold"),
#     plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
#     text = element_text(size=14)
#   )
# p1
# 
# 
# 
# #####
# 
# # null dataset
# 
# nonCDS <- subset(cv, !cv$H_ANN %like% "A")


### What is the potential for using mouse to understand the mechanistic function of pathogenic variants?

### are gwas snvs associated with disease more or less conserved than expected in mouse?




# read paper about human-mouse alignment methods
# read paper about Ensembl Regulatory Build


# p1.2 <- ggplot() +
#   geom_bar(data=plot_dt2, aes(x=H_annotation, y=M_ann_frac, fill=M_annotation), colour = "black", stat="identity", position=position_dodge()) +
#   xlab("Human genomic annotation") +
#   ylab("Mouse alignment composition (fraction)") +
#   ggtitle("B") +
#   labs(fill = "Mouse genomic\nannotation") +
#   coord_flip() +
#   scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),
#                      limits = c(0, 1)) +
#   scale_fill_brewer(palette="Set3") +
#   # scale_fill_brewer(palette="Paired") +
#   theme_bw() +
#   theme(
#     # legend.position = "none",
#     # legend.title = element_blank(),
#     legend.title = element_text(hjust = 0.5),
#     # legend.key.size = unit(2, 'lines'),
#     # legend.justification=c(1,0),
#     # legend.position=c(0.95, 0.05),
#     # legend.box.background = element_rect(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(size = 26, face = "bold"),
#     plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
#     text = element_text(size=14)
#   )
# p1.2
