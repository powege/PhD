rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)


cv <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_mouse_mapping.csv")
cv_z <- fread("~/Dropbox/PhD/Data/Interspecific_SNV_mapping/ClinVar_z_scores.csv")


colnames(cv) <- c("H_CHR", "H_POS", "H_REF", "H_Kmer", "H_ANN", "M_CHR", "M_POS", "M_REF", "M_Kmer", "M_ANN")
cv_z$H_annotation <- mapvalues(cv_z$H_annotation, from=c("B","C","D","E","F","G","H","I","J"), 
                                 to=c("Exon - UTR",
                                      "Exon - non-coding",
                                      "Promoter",
                                      "Enhancer",
                                      "Open chromatin",
                                      "TF binding",
                                      "Promoter flanking",
                                      "Intron",
                                      "Unannotated"))
p1_dt <- cv_z[,c("H_annotation", "F_align_z", "F_conserved_z", "F_conserved_K_z")]
p1_dt <- melt(p1_dt, id.vars = "H_annotation")
p1_dt$variable <- as.character(p1_dt$variable)
p1_dt$variable[p1_dt$variable == "F_align_z"] <- "Align to mouse"
p1_dt$variable[p1_dt$variable == "F_conserved_z"] <- "Conserved in\nmouse"
p1_dt$variable[p1_dt$variable == "F_conserved_K_z"] <- "Conserved 3-mer\nin mouse"


# mouse annotation composition:
annotation <- c("B", "C", "D", "E", "F", "G", "H", "I", "J")
H_annotation <- list()
M_annotation <- list()
H_ann_total <- list()
M_ann_total <- list()
cv_align <- cv[!is.na(cv$M_POS),]
cv_align <- subset(cv_align, !cv_align$H_ANN %like% "A")
for (i in 1:length(annotation)){
  sub <- subset(cv_align, cv_align$H_ANN %like% annotation[i])
  H_ann_total[[i]] <- rep(nrow(sub), 9)
  tmp_vec <- rep(NA, 9)
  for (j in 1:length(annotation)){
    tmp_vec[j] <- nrow(sub[sub$M_ANN %like% annotation[j],])
    # print(c(i,j))
  }
  M_ann_total[[i]] <- tmp_vec
  H_annotation[[i]] <- rep(annotation[i], 9)
  M_annotation[[i]] <- c("B", "C", "D", "E", "F", "G", "H", "I", "J")
}  

p2_dt <- data.table(
  H_annotation = unlist(H_annotation),
  H_ann_total = unlist(H_ann_total),
  M_annotation = unlist(M_annotation),
  M_ann_total = unlist(M_ann_total)
)
p2_dt$M_ann_frac <- p2_dt$M_ann_total/p2_dt$H_ann_total
p2_dt$H_annotation <- mapvalues(p2_dt$H_annotation, from=c("B","C","D","E","F","G","H","I","J"),
                                   to=c(
                                        "Exon - UTR",
                                        "Exon - non-coding",
                                        "Promoter",
                                        "Enhancer",
                                        "Open chromatin",
                                        "TF binding",
                                        "Promoter flanking",
                                        "Intron",
                                        "Unannotated"))
p2_dt$M_annotation <- mapvalues(p2_dt$M_annotation, from=c("B","C","D","E","F","G","H","I","J"),
                                   to=c(
                                        "Exon - UTR",
                                        "Exon - non-coding",
                                        "Promoter",
                                        "Enhancer",
                                        "Open chromatin",
                                        "TF binding",
                                        "Promoter flanking",
                                        "Intron",
                                        "Unannotated"))

# set factor order by total alignment
order_H <- c("Unannotated", "TF binding", "Intron", "Open chromatin", "Promoter flanking", "Enhancer",
             "Exon - non-coding", "Promoter", "Exon - UTR")
order_M <- rev(order_H)
p2_dt$H_annotation  <- factor(p2_dt$H_annotation, levels = as.character(order_H))
p2_dt$M_annotation  <- factor(p2_dt$M_annotation, levels = as.character(order_M))
p1_dt$H_annotation  <- factor(p1_dt$H_annotation, levels = as.character(order_H))

### PLOT

p1 <- ggplot(p1_dt, aes(x=H_annotation, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  # scale_y_continuous(trans = log2_trans()) +
  xlab("Human genomic annotation") +
  ylab("N pathogenic SNVs (z-score)") +
  ggtitle("A") +
  coord_flip() +
  theme_bw() +
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    legend.justification=c(1,0),
    legend.position=c(0.95, 0.65),
    # legend.box.background = element_rect(colour = "black"),
    # legend.spacing.x = unit(1.0, 'cm'),
    # legend.text = element_text(margin = margin(t = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) +
  guides(fill=guide_legend(
    keywidth=0.4,
    keyheight=0.4,
    default.unit="inch")
  )
p1

p2 <- ggplot() +
  geom_bar(data=p2_dt, aes(x=H_annotation, y=M_ann_frac, fill=M_annotation), colour = "black", stat="identity") +
  xlab("Human genomic annotation") +
  ylab("Mouse alignment composition (fraction)") +
  ggtitle("B") +
  labs(fill = "Mouse genomic\nannotation") +
  coord_flip() +
  # scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2.0, 2.5),
  #                    limits = c(0, 2.5)) +
  scale_fill_brewer(palette="Set3") +
  # scale_fill_brewer(palette="Paired") +
  theme_bw() +
  theme(
    # legend.position = "none",
    # legend.title = element_blank(),
    legend.title = element_text(hjust = 0.5),
    # legend.key.size = unit(2, 'lines'),
    # legend.justification=c(1,0),
    # legend.position=c(0.95, 0.05),
    # legend.box.background = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 26, face = "bold"),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  ) 
p2

pout <- grid.arrange(p1, p2, nrow = 1, widths = c(6, 8))
ggsave("~/Dropbox/PhD/Data/Figures/Figure_tmp4.jpg", plot = pout, height = 6, width = 14)




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
