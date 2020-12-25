### SCRIPT that:
    # Tests correlation in methylation intensity between ES and NPC. 
    # Plots fraction of CGs by methylation (%)
    # Plots substitution probabilities for CpGs across methylation groups (all, low, medium, high)

### TO DO
# ancestral mutant dataset
# remove QCed loci (GERP, coverage, sm)
# C > T mutations

rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)


meth <- fread("gunzip -cq ~/Dropbox/PhD/Data/Ensembl/Methylation/formatted/m.mus_GRC38_Ensembl_CG_meth_NPC_ES.bed.gz")
chr19 <- fread("~/Dropbox/PhD/Data/MGP/Variants/vcf_QCed_VEP/MGP_v5_allSTRAIN_snps_QCed_VEP_v94_chr19.vcf", 
              select = c(1:2),
              col.names = c("chromosome", "start"))
Pmu <- fread("~/Dropbox/PhD/Data/Mu_rates/MGP_v5_allSTRAIN_kmer_PSNV_specific.csv")

meth <- meth[coverage > 0]

all19 <- meth[chromosome == 19]
low19 <- meth[percent_methylated < 20 & chromosome == 19]
medium19 <- meth[percent_methylated >= 20 & percent_methylated <= 60 & chromosome == 19]
high19 <- meth[percent_methylated > 60 & chromosome == 19]

low19_mu <- low19[chr19, on = c("chromosome", "start")]
low19_mu <- low19_mu[complete.cases(low19_mu)]
low19_Pmu <- nrow(low19_mu) / nrow(low19)

medium19_mu <- medium19[chr19, on = c("chromosome", "start")]
medium19_mu <- medium19_mu[complete.cases(medium19_mu)]
medium19_Pmu <- nrow(medium19_mu) / nrow(medium19)

high19_mu <- high19[chr19, on = c("chromosome", "start")]
high19_mu <- high19_mu[complete.cases(high19_mu)]
high19_Pmu <- nrow(high19_mu) / nrow(high19)

CG_Pmu <- Pmu$k1CG_mu_rate[Pmu$k1CG == "C (CG)" & Pmu$to == "T"][1]
nonCG_Pmu <- Pmu$k1CG_mu_rate[Pmu$k1CG == "C (nonCG)" & Pmu$to == "T"][1]


#####


p1 <- ggplot(meth, aes(x=percent_methylated)) + 
  geom_histogram(binwidth=10, colour = "black", fill = "green") +
  xlab("Methylation (%)") +
  ylab("Count") +
  geom_vline(xintercept=20, linetype="dashed", color = "blue", size=1) +
  geom_vline(xintercept=60, linetype="dashed", color = "blue", size=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p1
ggsave(filename = "~/Dropbox/tmp_meth_hist.jpg", plot = p1, height = 4, width = 5)

bar1_dt <- data.table(group = c("CG", "non-CG"),
                      pMU = c(CG_Pmu, nonCG_Pmu))
p2 <- ggplot(data=bar1_dt, aes(x=group, y=pMU, fill=group)) +
  geom_bar(stat="identity", colour="black") +
  xlab("") +
  ylab("P(substitution)") +
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p2
ggsave(filename = "~/Dropbox/tmp_meth_bar1.jpg", plot = p2, height = 5, width = 5)

bar2_dt <- data.table(group = c("High", "Medium", "Low"),
                      pMU = c(high19_Pmu, medium19_Pmu, low19_Pmu))
p3 <- ggplot(data=bar2_dt, aes(x=group, y=pMU, fill=group)) +
  geom_bar(stat="identity", colour="black") +
  xlab("") +
  ylab("P(substitution)") +
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.1, 0.15, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm"),
    text = element_text(size=14)
  )
p3
ggsave(filename = "~/Dropbox/tmp_meth_bar2.jpg", plot = p3, height = 3, width = 3)






