rm(list=ls())
graphics.off()

library(ggplot2)
library(scales)

list.files("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/", pattern = "FRACTION")

PC_exon <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_PC_exon_POS_short.csv", header = F)
enhancer <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_enhancer_POS_short.csv", header = F)
TF <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_TF_binding_POS_short.csv", header = F)
NC_exon <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_NC_exon_POS_short.csv", header = F)
promoter <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_promoter_POS_short.csv", header = F)
UTR <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_UTR_POS_short.csv", header = F)
OC <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_open_chromatin_POS_short.csv", header = F)
intron <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_intron_POS_short.csv", header = F)
unannotated <- read.csv("~/Dropbox/PhD/Data/NC_constraint/Ensembl_annotation_POS/FRACTION_unannotated_POS_short.csv", header = F)

colnames(PC_exon) <- c("percentile", "fraction")
colnames(enhancer) <- c("percentile", "fraction")
colnames(TF) <- c("percentile", "fraction")
colnames(NC_exon) <- c("percentile", "fraction")
colnames(promoter) <- c("percentile", "fraction")
colnames(UTR) <- c("percentile", "fraction")
colnames(OC) <- c("percentile", "fraction")
colnames(intron) <- c("percentile", "fraction")
colnames(unannotated) <- c("percentile", "fraction")

PC_exon$annotation <- "Exon, CDS"
enhancer$annotation <- "Enhancer"
TF$annotation <- "TF binding"
NC_exon$annotation <- "Exon, ncRNA"
promoter$annotation <- "Promoter"
UTR$annotation <- "Exon, UTR"
OC$annotation <- "Open chromatin"
intron$annotation <- "Intron"
unannotated$annotation <- "Unannotated"

# calculate fold change vs the 100th percentile
PC_exon$fold_change <- PC_exon$fraction/PC_exon$fraction[PC_exon$percentile == 100]
enhancer$fold_change <- enhancer$fraction/enhancer$fraction[enhancer$percentile == 100]
TF$fold_change <- TF$fraction/TF$fraction[TF$percentile == 100]
NC_exon$fold_change <- NC_exon$fraction/NC_exon$fraction[NC_exon$percentile == 100]
promoter$fold_change <- promoter$fraction/promoter$fraction[promoter$percentile == 100]
UTR$fold_change <- UTR$fraction/UTR$fraction[UTR$percentile == 100]
OC$fold_change <- OC$fraction/OC$fraction[OC$percentile == 100]
intron$fold_change <- intron$fraction/intron$fraction[intron$percentile == 100]
unannotated$fold_change <- unannotated$fraction/unannotated$fraction[unannotated$percentile == 100]

df <- rbind(NC_exon, promoter, UTR, OC, enhancer, PC_exon, TF, intron, unannotated)

# df$Outlier <- 0
# df$Outlier[df$Percentile == 1] <- 1
# df$Outlier <- as.factor(df$Outlier)

p1 <- ggplot(df, aes(x=percentile, y=fold_change, color=annotation)) + 
  # geom_point() +
  stat_smooth(se = F) +
  # scale_color_manual(values=c('grey', 'red')) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
scale_y_continuous(trans = log2_trans(),
                   # breaks = c(1, 2, 4, 8, 16)) +
                    breaks = trans_breaks("log2", function(x) 2^x)) +
scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14)
  )
p1
ggsave(filename = "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_annotation_fold_change_by_constraint.jpeg",
       plot = p1,
       height = 4, 
       width = 6)

###########

gobbler <- function(promoter){
  percentile_group <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                        "51-75", "76-90", "91-95", "96-98", "99", "100")
  mean_fold_change <- rep(NA, 12)
  mean_fold_change[1] <- mean(promoter$fold_change[promoter$percentile %in% c(1)])
  mean_fold_change[2] <- mean(promoter$fold_change[promoter$percentile %in% c(2)])
  mean_fold_change[3] <- mean(promoter$fold_change[promoter$percentile %in% c(3:5)])
  mean_fold_change[4] <- mean(promoter$fold_change[promoter$percentile %in% c(6:10)])
  mean_fold_change[5] <- mean(promoter$fold_change[promoter$percentile %in% c(11:25)])
  mean_fold_change[6] <- mean(promoter$fold_change[promoter$percentile %in% c(26:50)])
  mean_fold_change[7] <- mean(promoter$fold_change[promoter$percentile %in% c(51:75)])
  mean_fold_change[8] <- mean(promoter$fold_change[promoter$percentile %in% c(76:90)])
  mean_fold_change[9] <- mean(promoter$fold_change[promoter$percentile %in% c(91:95)])
  mean_fold_change[10] <- mean(promoter$fold_change[promoter$percentile %in% c(96:98)])
  mean_fold_change[11] <- mean(promoter$fold_change[promoter$percentile %in% c(99)])
  mean_fold_change[12] <- mean(promoter$fold_change[promoter$percentile %in% c(100)])
  out <- data.frame(percentile_group = percentile_group,
                    rank = c(1:12),
                    mean_fold_change = mean_fold_change,
                    annotation = rep(promoter$annotation[1], 12))
  return(out)
}

PC_exon2 <- gobbler(PC_exon)
enhancer2 <- gobbler(enhancer)
TF2 <- gobbler(TF)
NC_exon2 <- gobbler(NC_exon)
promoter2 <- gobbler(promoter)
UTR2 <- gobbler(UTR)
OC2 <- gobbler(OC)
intron2 <- gobbler(intron)
unannotated2 <- gobbler(unannotated)

df2 <- rbind(NC_exon2, promoter2, UTR2, OC2, enhancer2, PC_exon2, TF2, intron2, unannotated2)
levels(df2$percentile_group) <- c("1", "2", "3-5", "6-10", "11-25", "26-50",
                                  "51-75", "76-90", "91-95", "96-98", "99", "100")

p2 <- ggplot(df2, aes(x=rank, y=mean_fold_change, color=annotation)) + 
  # geom_point() +
  geom_line(size=1.5) +
  xlab("Constraint (percentile rank)") +
  ylab("Genomic territory\n(fold change versus 100th percentile)") +
  scale_y_continuous(trans = log2_trans(),
                     # breaks = trans_breaks("log2", function(x) 2^x)) +
                     breaks = c(0.25, 0.5, 1, 2, 5, 10, 20, 40),
                     limits = c(0.25, 40)) +
  scale_x_continuous(breaks=c(1:12),
                     labels=c("1", "2", "3-5", "6-10", "11-25", "26-50",
                              "51-75", "76-90", "91-95", "96-98", "99", "100")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p2
ggsave(filename = "~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Figure_annotation_fold_change_by_constraint_2.jpeg",
       plot = p2,
       height = 5, 
       width = 7)

#######

# ggplot(NC_exon, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(OC, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(promoter, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(UTR, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(PC_exon, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(enhancer, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(TF, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 
# 
# ggplot(intron, aes(x=percentile, y=fold_change)) + 
#   geom_point() +
#   stat_smooth() +
#   # scale_color_manual(values=c('grey', 'red')) +
#   xlab("Constraint percentile") +
#   ylab("Genomic territory (fraction)") 


#######
