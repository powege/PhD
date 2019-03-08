rm(list = ls())
graphics.off()

library(data.table)
# library(dplyr)

# set variables
contraint_variables_path <- "~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_variables_by_window.csv"
output_file_path <- "~/Dropbox/PhD/Data/NC_constraint/MGP_constraint_by_window.csv"

# import
dt <- fread(contraint_variables_path)

# QC
# summary(dt)
dt <- dt[complete.cases(dt),]

# Calculate constraint
mod <- lm(dt$n_SNV ~ dt$p_SNV_given_kmers + dt$Read_depth)
# summary(mod)
dt$Constraint_score <- studres(mod)

# percentile rank
percentile_rank <- function(x) ceiling((rank(x)/length(x))*100) 
dt$Constraint_percentile <- percentile_rank(dt$Constraint_score)


# Output
output <- dt[, c("CHR", "POS_from", "POS_to", "Constraint_score", "Constraint_percentile")]
fwrite(output, output_file_path)



#####

sub.int <- sample(nrow(dt), 100000, replace = F)
sub <- dt[sub.int,]
plot(sub$n_SNV ~ sub$p_SNV_given_kmers)
abline(lm(sub$n_SNV ~ sub$p_SNV_given_kmers), col = "red")
abline(lm(sub$n_SNV ~ sub$p_SNV_given_kmers + sub$Read_depth), col = "blue")
# plot(lm(sub$n_SNV ~ sub$p_SNV_given_kmers + sub$Read_depth))


plot(sub$n_SNV ~ sub$Read_depth)
abline(lm(sub$n_SNV ~ sub$Read_depth), col = "red")
