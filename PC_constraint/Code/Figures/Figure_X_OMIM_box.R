### FIGURE 3 -- 

rm(list = ls())
graphics.off()

library(data.table)
library(ggplot2)
library(plotROC)
library(pROC)
library("cowplot")
library(Hmisc)
library(ggpubr)

### IMPORT 

omim <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Disease/OMIM_Petrovski2013_Dataset_S1.csv")
newC <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/HM_constriant_orth.csv")
RVIS_esp <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/RVIS_Petrovski2013_ESP.csv")
Lek <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Constraint_scores/H_CS_fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt")


### FORMAT

### # convert "" to NA
omim[omim==""] <- NA

# colnames
colnames(RVIS_esp) <- c("H_external_gene_name", "RVIS_ESP", "RVIS_ESP_percentile")
Lek <- Lek[,c("gene", "mis_z", "lof_z", "pLI")]
colnames(Lek) <- c("H_external_gene_name", "mis_z_Lek", "lof_z_Lek", "pLI_Lek")

# percentiles
# RVIS_esp$RVIS_ESP_percentile <- RVIS_esp$RVIS_ESP_percentile/100
# percentile <- ecdf(Lek$mis_z_Lek)
# Lek$mis_z_percentile_Lek <- percentile(Lek$mis_z_Lek)
# percentile <- ecdf(Lek$pLI_Lek)
# Lek$pLI_percentile_Lek <- percentile(Lek$pLI_Lek)
# percentile <- ecdf(newC$H_fun_Z_0.001)
# newC$H_fun_Z_percentile <- percentile(newC$H_fun_Z_0.001)
# percentile <- ecdf(newC$M_fun_Z_0.0001)
# newC$M_fun_Z_percentile <- percentile(newC$M_fun_Z_0.0001)

# df <- merge(Lek, RVIS_esp, all = T)
# df <- merge(df, newC, all = T)
df <- newC[,c("H_external_gene_name", "M_external_gene_name",
              "H_fun_Z_0.001", "M_fun_Z_0.0001")]
df <- df[!is.na(df$H_external_gene_name),]
df <- df[df$H_external_gene_name != "",]

# create dataframes gene names and OMIM category, inlcuding all genes
recesive <- data.frame(H_external_gene_name = na.omit(unique(omim$`OMIM Recessive`)),
                       cat = rep("Recessive", length(na.omit(unique(omim$`OMIM Recessive`)))))
                       
haploinsufficient <- data.frame(H_external_gene_name = na.omit(unique(omim$`OMIM Haploinsufficiency`)),
                       cat = rep("Haploinsufficient", length(na.omit(unique(omim$`OMIM Haploinsufficiency`)))))

dom.neg <- data.frame(H_external_gene_name = na.omit(unique(omim$`OMIM Dominant Negative`)),
                       cat = rep("Dominant-\nnegative", length(na.omit(unique(omim$`OMIM Dominant Negative`)))))

denovo <- data.frame(H_external_gene_name = na.omit(unique(omim$`OMIM de novo`)),
                       cat = rep("De novo", length(na.omit(unique(omim$`OMIM de novo`)))))

omim.all <- data.frame(H_external_gene_name = na.omit(unique(omim$`OMIM disease genes`)),
                       cat = rep("All OMIM", length(na.omit(unique(omim$`OMIM disease genes`)))))

non.OMIM <- df$H_external_gene_name
non.OMIM <- non.OMIM[!non.OMIM %in% omim.all$H_external_gene_name]
non.OMIM <- data.frame(H_external_gene_name = non.OMIM,
                       cat = rep("non-OMIM", length(non.OMIM)))

df.cat <- rbind(recesive, haploinsufficient, dom.neg, denovo, omim.all, non.OMIM)
df <- df[df.cat, on = "H_external_gene_name"]
df <- df[!duplicated(df),]
df <- df[complete.cases(df),]
table(df$cat)

# percentiles
percentile <- ecdf(df$M_fun_Z_0.0001)
df$M_fun_Z_percentile <- percentile(df$M_fun_Z_0.0001)

M_n.hap <- length(df$M_external_gene_name[df$cat == "Haploinsufficient"])
M_n.rec <- length(df$M_external_gene_name[df$cat == "Recessive"])
M_n.DN <- length(df$M_external_gene_name[df$cat == "Dominant-\nnegative"])
M_n.dvo <- length(df$M_external_gene_name[df$cat == "De novo"])
M_n.OMIM <- length(df$M_external_gene_name[df$cat == "All OMIM"])
M_n.all <- length(df$M_external_gene_name[df$cat == "non-OMIM"])

df$cat <- as.character(df$cat)
df$cat[df$cat == "Haploinsufficient"] <- paste0("Haploinsufficient\n(n=", M_n.hap, ")")
df$cat[df$cat == "Recessive"] <- paste0("Recessive\n(n=", M_n.rec, ")")
df$cat[df$cat == "Dominant-\nnegative"] <- paste0("Dominant-\nnegative\n(n=", M_n.DN, ")")
df$cat[df$cat == "De novo"] <- paste0("De novo\n(n=", M_n.dvo, ")")
df$cat[df$cat == "All OMIM"] <- paste0("All OMIM\n(n=", M_n.OMIM, ")")
df$cat[df$cat == "non-OMIM"] <- paste0("Non-OMIM\n(n=", M_n.all, ")")
table(df$cat)


df$cat <- reorder(df$cat, df$cat, FUN = length, order = T)
table(df$cat)
comparisons <- list( c("Non-OMIM\n(n=14176)", "Haploinsufficient\n(n=156)"), 
                       c("Non-OMIM\n(n=14176)", "Dominant-\nnegative\n(n=328)"),
                       c("Non-OMIM\n(n=14176)", "De novo\n(n=416)"),
                     c("Non-OMIM\n(n=14176)", "Recessive\n(n=743)"),
                     c("Non-OMIM\n(n=14176)", "All OMIM\n(n=2094)"))

box <- ggplot(df, aes(x=cat, y=M_fun_Z_percentile, fill = cat)) +
  geom_boxplot() +
  stat_compare_means(
    aes(group = cat),
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    label.y = c(1.05, 1.15, 1.25, 1.35, 1.45),
    paired = F
    # size = 3,
  ) +
  ylab("Mouse orthologue funZ percentile") +
  xlab("") +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.key.size = unit(1.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1),"cm")) 
box

ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_OMIM_box.jpg", plot = box, width = 9, height = 7)


############
# 
# # Set gene lists
# recesive <- na.omit(unique(omim$`OMIM Recessive`))
# df$OMIM.recessive <- 0
# df$OMIM.recessive[df$H_external_gene_name %in% recesive] <- 1
# 
# haploinsufficient <- na.omit(unique(omim$`OMIM Haploinsufficiency`))
# df$OMIM.hap <- 0
# df$OMIM.hap[df$H_external_gene_name %in% haploinsufficient] <- 1
# 
# dom.neg <- na.omit(unique(omim$`OMIM Dominant Negative`))
# df$OMIM.domneg <- 0
# df$OMIM.domneg[df$H_external_gene_name %in% dom.neg] <- 1
# 
# denovo <- na.omit(unique(omim$`OMIM de novo`))
# df$OMIM.denovo <- 0
# df$OMIM.denovo[df$H_external_gene_name %in% denovo] <- 1
# 
# omim.all <- na.omit(unique(omim$`OMIM disease genes`))
# df$OMIM.all <- 0
# df$OMIM.all[df$H_external_gene_name %in% omim.all] <- 1
# 
# omim.all <- na.omit(unique(omim$`OMIM disease genes`))
# df$OMIM.all <- 0
# df$OMIM.all[df$H_external_gene_name %in% omim.all] <- 1
# 
# # subset one to one orthologs
# # df <- subset(df, df$orthology_type == "ortholog_one2one")
# 
# ### OMIM Boxplot
# x <-c("M_fun_Z_percentile_0.0001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all")
# boxplot.OMIM <- function(x, df){
#   
#   box.df <- df[ ,x, with = F]
#   box.df <- box.df[which(!is.na(box.df$M_fun_Z_percentile_0.0001)),]
#   n.recessive <- length(box.df$OMIM.recessive[box.df$OMIM.recessive == 1])
#   n.hap <- length(box.df$OMIM.hap[box.df$OMIM.hap == 1])
#   n.domneg <- length(box.df$OMIM.domneg[box.df$OMIM.domneg == 1])
#   n.denovo <- length(box.df$OMIM.denovo[box.df$OMIM.denovo == 1])
#   n.all <- length(box.df$OMIM.all[box.df$OMIM.all == 1])
#   colnames(box.df) <- c(x[1], 
#                         paste0("Recessive\n(n=", n.recessive, ")"),
#                         paste0("Haploinsufficient\n(n=", n.hap, ")"),
#                         paste0("Dominant-\nnegative\n(n=", n.domneg, ")"),
#                         paste0("De novo\n(n=", n.denovo, ")"),
#                         paste0("All\n(n=", n.all, ")"))
#   box.df <- melt(box.df, id.vars = x[1])
#   box.df$value[box.df$value == 1] <- "In OMIM list"
#   box.df$value[box.df$value == 0] <- "Not in OMIM list"
#   ind <- which(box.df[,1] > -2 & box.df[,1] < 2)
#   box.df <- box.df[ind,]
#   box.df <- as.data.frame(box.df)
#   box.df$value <- as.factor(box.df$value)
#   out.box <- ggplot(box.df, aes_string(x=box.df[,2], y=box.df[,1], fill=box.df[,3])) +
#     # geom_violin() +
#     geom_boxplot() +
#     # geom_boxplot(width=0.1, position=position_dodge(0.9)) +
#     stat_compare_means(
#       aes(group = box.df[,3]),
#       method = "wilcox.test", 
#       label = "p.signif",
#       paired = F,
#       label.y = c(1),
#       size = 7
#     ) +
#     ylab("Mouse funZ percentile") +
#     xlab("") +
#     scale_fill_manual(values=c("plum3", "white")) +
#     theme_classic() +
#     theme(legend.title=element_blank(),
#           axis.text = element_text(size = 14),
#           axis.title = element_text(size = 14),
#           legend.text=element_text(size=14),
#           legend.position="top",
#           plot.margin=unit(c(1,1,1,1),"cm")) 
#   out.box
#   # medians <- box.df[, .(median(value)), by = c("IMPC.Viability.2", "variable")][order(IMPC.Viability.2, variable)]
#   return(out.box)
# }
# 
# plotA <- boxplot.OMIM(c("M_fun_Z_percentile_0.0001", "OMIM.recessive", "OMIM.hap", "OMIM.domneg", "OMIM.denovo", "OMIM.all"), df)
# plotA
# ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Figures/Figure_3A.jpg", plot = plotA, width = 8, height = 6)
