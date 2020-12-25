rm(list = ls())
graphics.off()

library(data.table)

source("~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R")

ann <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz")

ann <- ann[,c("chromosome",
              "start",
              "end",
              "category",
              "category_id",
              "regulatory_feature_stable_id",
              "transcript_id",
              "percentage_Nmask",
              "percentage_sm",
              "strand")]

test <- ann[category == "Exon - UTR" | category == "Promoter"]

train_pop <- ann[category == "Intron - distal" | category == "Unannotated"]
train_pop <- collapse.overlap(train_pop[,1:3])


test_length <- (test$end + 1) - test$start
hist(test_length)


### STACK OVERFLOW

# how to efficienntly sample sequences from a bed file?

# I need to sample sequences of differing length from a bed file (a file that privides the start and end coordinnates of a sequence, and a category). 
# For example, given the bed file: 
bed <- data.table(category = c("A", "A", "A", "A", "B", "B"),
                  start = c(1, 100, 300, 410, 1, 810),
                  end = c(80, 220, 400, 700,  400, 900))

# And a vector of sequence lengths
seq_lengths <- c(7, 5, 3, 6, 6, 10)

# How can i randomly sample the same number and length of sequences from seq_length
# from within the bed file coordinates? The output would be in bed file format, something like:
sample <- data.table(category = c("A", "A", "A", "A", "B", "B"),
                     start = )

# The dataset I am applying this to is very large, and so performance is important.


# For example, i can do this by collapsing the sequences into a large vector and sampling.

# x = population sequence start
# y = population sequence end
# p = sample sequence length 
extractRandWindow <- function(x, y, p){
  firstIndex = sample(seq(((y+1)-x) - p + 1), 1)
  c(x:y)[firstIndex:(firstIndex + p -1)]
}
extractRandWindow(x=1, y=100, p=10)

y = 100
x = 1
p = 10

