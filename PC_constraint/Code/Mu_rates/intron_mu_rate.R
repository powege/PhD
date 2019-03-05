### SCRIPT that calculates the intron mutation rate for each Ensembl transcript

rm(list = ls())
graphics.off()

library(data.table)
library(plyr)

### FUNCTION that sums intron length from exon start and end POS
intron_len <- function(sub){
  dif <- sub$exon_chrom_end - sub$exon_chrom_start
  len <- sum(dif) + length(dif)
  out <- c(sub$ensembl_transcript_id[1], len)
  names(out) <- c("ensembl_transcript_id", "intron_length")
  return(out)
}

# HUMAN
dt <- fread("../../Data/Ensembl/H_exon_POS.csv")
output <- ddply(dt, "ensembl_transcript_id", intron_len)
fwrite(output, "../../Data/Ensembl/H_intron_length.csv")

# MOUSE
dt <- fread("../../Data/Ensembl/M_exon_POS.csv")
output <- ddply(dt, "ensembl_transcript_id", intron_len)
fwrite(output, "../../Data/Ensembl/M_intron_length.csv")


