rm(list=ls())
graphics.off()

library(data.table)
library(ggplot2)

### IMPORT 

# chr1
# five <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr1_500_100.csv")
# seven <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr1_700_100.csv")
# nine <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr1_900_100.csv")
# eleven <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr1_1100_100.csv")

# all chr 
five.in <- list()
for (chr in 1:19){
  five.in[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", chr, "_500_100.csv"))
}
five <- rbind.fill(five.in)
rm(five.in)
seven.in <- list()
for (chr in 1:19){
  seven.in[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", chr, "_700_100.csv"))
}
seven <- rbind.fill(seven.in)
rm(seven.in)
nine.in <- list()
for (chr in 1:19){
  nine.in[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", chr, "_900_100.csv"))
}
nine <- rbind.fill(nine.in)
rm(nine.in)
eleven.in <- list()
for (chr in 1:19){
  eleven.in[[chr]] <- fread(paste0("~/Dropbox/PhD/Data/NC_constraint/Constraint/MGP_constraint_variables_by_window_chr", chr, "_1100_100.csv"))
}
eleven <- rbind.fill(eleven.in)
rm(eleven.in)

### EXAMINE DATA

# histogram repeats
hist(five$Repeats, breaks = 100)
hist(seven$Repeats, breaks = 100)
hist(nine$Repeats, breaks = 100)
hist(eleven$Repeats, breaks = 100)

# histogram read depth
hist(five$Read_depth, breaks = 100)
hist(seven$Read_depth, breaks = 100)
hist(nine$Read_depth, breaks = 100)
hist(eleven$Read_depth, breaks = 100)

# histogram p_SNV
hist(five$p_SNV_given_kmers, breaks = 100)
hist(seven$p_SNV_given_kmers, breaks = 100)
hist(nine$p_SNV_given_kmers, breaks = 100)
hist(eleven$p_SNV_given_kmers, breaks = 100)

### FORMAT

# make dataframes match
# seven <- seven[3:(nrow(seven)-2),]
# nine <- nine[2:(nrow(nine)-1),]

# remove rows with NA
five <- five[complete.cases(five),]
seven <- seven[complete.cases(seven),]
nine <- nine[complete.cases(nine),]
eleven <- eleven[complete.cases(eleven),]

# remove regions with poor coverage
five.cut <- five[five$Read_depth >= 0.5,]
seven.cut <- seven[seven$Read_depth >= 0.5,]
nine.cut <- nine[nine$Read_depth >= 0.5,]
eleven.cut <- eleven[eleven$Read_depth >= 0.5,]

# remove highly repetative regions
five.cut <- five.cut[five.cut$Repeats < 1,]
seven.cut <- seven.cut[seven.cut$Repeats < 1,]
nine.cut <- nine.cut[nine.cut$Repeats < 1,]
eleven.cut <- eleven.cut[eleven.cut$Repeats < 1,]

hist(five.cut$Repeats, breaks = 500)
hist(seven.cut$Repeats, breaks = 700)
hist(nine.cut$Repeats, breaks = 900)
hist(eleven.cut$Repeats, breaks = 1100)

hist(five.cut$Read_depth, breaks = 500)
hist(seven.cut$Read_depth, breaks = 700)
hist(nine.cut$Read_depth, breaks = 900)
hist(eleven.cut$Read_depth, breaks = 1100)

hist(five.cut$n_SNV, breaks = 500)
hist(seven.cut$n_SNV, breaks = 700)
hist(nine.cut$n_SNV, breaks = 900)
hist(eleven.cut$n_SNV, breaks = 1100)

