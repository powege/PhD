rm(list=ls())
graphics.off()

library(data.table)

### IMPORT
mk3 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_3mer_mu_rate.table")
hk3 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_3mer_mu_rate.table")
mk5 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_5mer_mu_rate.table")
hk5 <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table")

mk3p <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/M_k3mer_canPC_Pmu.csv")
hk3p <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k3mer_canPC_Pmu.csv")
hk5p <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Sequence_Pmu/H_k5mer_canPC_Pmu.csv")



### SUMMARY 
summary(mk3)
summary(hk3)
summary(mk5)
summary(hk5)

### CORRELATION
plot(hk3$k3_mu_rr, mk3$k3_mu_rr)
cor.test(hk3$k3_mu_rr, mk3$k3_mu_rr)
plot(hk5$k5_mu_rr, mk5$k5_mu_rr)
cor.test(hk5$k5_mu_rr, mk5$k5_mu_rr)
plot(hk3p$p_syn, hk5p$p_syn)
cor.test(hk3p$p_syn, hk5p$p_syn)
plot(hk3p$p_mis, hk5p$p_mis)
cor.test(hk3p$p_mis, hk5p$p_mis)


