rm(list=ls())
graphics.off()

library(data.table)
library(plyr)

M_nvar <- rep(NA, 20)
chr <- c(1:19)
for (i in 1:length(chr)){
  dt <- fread(paste0("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/M_MGP_QCed_VEP_all_chr", chr[i], ".txt"))
  M_nvar[i] <- nrow(dt)
}
M_ntot <- sum(M_nvar)

H_nvar_0.0001 <- rep(NA, 23)
H_nvar_0.0005 <- rep(NA, 23)
H_nvar_0.001 <- rep(NA, 23)
chr <- c(1:22)
for (i in 1:length(chr)){
  dt <- fread(paste0("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/VEP/Formatted_output/H_1000GP_QCed_VEP_all_chr", chr[i], ".txt"), fill = T)
  dt.0005 <- dt[dt$V19 > 0.0005,]
  H_nvar_0.0005[i] <- nrow(dt.0005)
  dt.001 <- dt[dt$V19 > 0.001,]
  H_nvar_0.001[i] <- nrow(dt.001)
  dt.0001 <- dt[dt$V19 > 0.0001,]
  H_nvar_0.0001[i] <- nrow(dt.0001)
}
H_ntot_0.0005 <- sum(H_nvar_0.0005)
H_ntot_0.001 <- sum(H_nvar_0.001)
H_ntot_0.0001 <- sum(H_nvar_0.0001)

x <- head(dt.001$V2)
diff(x)
x <- diff(dt.0005$V2)
x <- x[x<1000]
hist(x, breaks = 100)

b <- x[x>600]

