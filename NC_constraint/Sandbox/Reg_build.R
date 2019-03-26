rm(list=ls())
graphics.off()

library(dplyr)
library(tidyr)
library(data.table)

x <-  fread("~/Dropbox/PhD/Data/Ensembl/Regulatory_build/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff")
colnames(x) <- c("chromosome", "source", "feature", "satrt", "end", "score", "strand", "frame", "attribute")
x <- 
 x %>% separate(attribute, c("ID", "bound_end","bound_start", "description", "feature_type"), sep = ";")
x$ID <- gsub("ID=", "", x$ID)
x$bound_end <- gsub("bound_end=", "", x$bound_end)
x$bound_start <- gsub("bound_start=", "", x$bound_start)
x$description <- gsub("description=", "", x$description)
x$feature_type <- gsub("feature_type=", "", x$feature_type)

table(x$feature_type)
table(x$description)

y <- fread("~/Dropbox/PhD/Data/GENCODE/gencode.vM20.annotation.gff3", fill = T)
colnames(y) <- c("chromosome", "source", "type", "satrt", "end", "score", "strand", "phase", "attributes")
head(y)
table(y$type)

test <- subset(y, y$type == "CDS")
head(test)
grep("gene_type=", y$attributes)





