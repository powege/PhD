rm(list = ls())
graphics.off()

library(data.table)
library(MASS)
# library(dplyr)

# set variables
contraint_variables_path <- "~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_variables_by_window.csv"
output_file_path <- "~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_by_window.csv"

# import
dt <- fread(contraint_variables_path)

# QC
# summary(dt)
dt <- dt[complete.cases(dt),]

# remove windows with read depth < 50%
len_dt_all <- nrow(dt)
dt1 <- subset(dt, dt$Read_depth >= 0.5)
repeat_windows <- subset(dt1, dt1$Repeats > 0.8)
dt1 <- subset(dt1, dt1$Repeats <= 0.8)

# Calculate constraint
mod <- lm(dt1$n_SNV ~ dt1$p_SNV_given_kmers
          # + dt1$Read_depth
          # + dt1$Repeats
          )
summary(mod)
dt1$Constraint_score <- studres(mod)

# percentile rank
percentile_rank <- function(x) ceiling((rank(x)/length(x))*100) 
dt1$Constraint_percentile <- percentile_rank(dt1$Constraint_score)

# Output
output <- dt1[, c("CHR", "POS_from", "POS_to", "Constraint_score", "Constraint_percentile")]
fwrite(output, output_file_path)



#####

# sub.int <- sample(nrow(dt), 100000, replace = F)
# sub <- dt[sub.int,]
# plot(sub$n_SNV ~ sub$p_SNV_given_kmers)
# abline(lm(sub$n_SNV ~ sub$p_SNV_given_kmers), col = "red")
# abline(lm(sub$n_SNV ~ sub$p_SNV_given_kmers + sub$Read_depth), col = "blue")
# # plot(lm(sub$n_SNV ~ sub$p_SNV_given_kmers + sub$Read_depth))
# 
# 
# plot(sub$n_SNV ~ sub$Read_depth)
# abline(lm(sub$n_SNV ~ sub$Read_depth), col = "red")


hist(repeat_windows$n_SNV)

hist(dt1$n_SNV, breaks = 100)
hist(dt1$Repeats)

per1 <- subset(dt1, dt1$Constraint_percentile == 1)
per2 <- subset(dt1, dt1$Constraint_percentile == 2)
per3 <- subset(dt1, dt1$Constraint_percentile == 3)
p1to3 <- subset(dt1, dt1$Constraint_percentile == 1 |
                  dt1$Constraint_percentile == 2 |
                  dt1$Constraint_percentile == 3)

hist(per1$n_SNV, breaks = 100)
hist(per2$n_SNV, breaks = 100)
hist(per3$n_SNV, breaks = 100)
hist(p1to3$n_SNV, breaks = 100)
hist(per1$Repeats)
hist(dt1$Read_depth)
hist(per1$Read_depth)
hist(per2$Read_depth)
hist(per3$Read_depth)

p1to5 <- dt1[dt1$Constraint_percentile %in% c(1:10),]
hist(p1to5$n_SNV, breaks = 100)
hist(p1to5$Repeats)

