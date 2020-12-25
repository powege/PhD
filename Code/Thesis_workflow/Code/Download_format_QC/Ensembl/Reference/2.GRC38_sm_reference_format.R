### INPUT 
# raw fasta files -- .fa.gz
# chromosome (as.character)

### OUTPUT
# N mask pos (cols: chromosome; start; end) -- .csv
# Soft mask pos (cols: chromosome; start; end) -- .csv.gz
# CG pos (cols: chromosome; start; end) -- .csv.gz
# bp counts (cols: chromosome; start; end; n_N_masked) -- .csv
# reference vector (cols: reference) -- .csv.gz

rm(list = ls())
graphics.off()

library(data.table)
library(stringr)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

# set args variables
in.file.fasta <- args[1]
chr <- as.character(args[2])
functions.file <- args[3]

out.file.N <- args[4]
out.file.sm <- args[5]
out.file.CG <- args[6]
out.file.bp <- args[7]
out.file.ref <- args[8]

# in.file.fasta <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/raw/Mus_musculus.GRCm38.dna_sm.chromosome.19.fa.gz"
# chr <- as.character(19)
# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Data_format_and_QC/Misc/FUNCTIONS.R"

# chr = "19"
# functions.file <- "/well/lindgren/George/Code/Workflows/FUNCTIONS.R"
# in.file.fasta <- "/well/lindgren/George/Data/Mapping_v1/raw/Mus_musculus.GRCm38.dna_sm.chromosome.19.fa.gz"
# out.file.N <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/m.mus_reference_Ensembl_GRC38_N_mask.bed.gz"
# out.file.sm <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/m.mus_reference_Ensembl_GRC38_soft_mask.bed.gz"
# out.file.CG <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/m.mus_reference_Ensembl_GRC38_CG.bed.gz"
# out.file.bp <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/m.mus_reference_Ensembl_GRC38_bp_counts.bed"
# out.file.ref <- "/well/lindgren/George/Data/Ensembl/Reference/Formatted/m.mus_reference_Ensembl_GRC38_chr19.csv.gz"

### FUNCTIONS

source(functions.file)

### IMPORT

fasta <- fread(paste0("gunzip -cq ", in.file.fasta), sep = "\t") # unnzip and import raw fasta

### FORMAT

meta <- names(fasta) # get header description
colnames(fasta) <- "V1" # set colname

meta.split <- strsplit(meta, ":") # identify satrt and end POS from header line
start.POS <- as.integer(meta.split[[1]][5])
end.POS <- as.integer(meta.split[[1]][6])

POS <- seq(from = start.POS, to = end.POS) # create vectors for POS and REF
REF <- paste(fasta$V1, collapse = '')
REF <- strsplit(REF, "")[[1]]

dt_ref <- data.table(POS = POS, 
                     REF = REF)

is.sequential <- function(x){ all(abs(diff(x)) == 1) }  # check sequence is continuous
if (dt_ref$POS[1]!=1) { stop("Start not equal to one", call.=FALSE) } 
if (!is.sequential(dt_ref$POS)) { stop("POS not sequential", call.=FALSE) } 

# N coordinates bed file
dt_N <- dt_ref[REF == "N" | REF == "n"] # subset N coordinates
dt_N$chromosome <- chr
N_out <- long.to.bed(dt_N[,c("chromosome", "POS")])

# sof masked coordinates bed file
dt_sm <- dt_ref[REF == "N" |
                  REF == "n" |
                  REF == "a" | 
                  REF == "c" | 
                  REF == "t" |
                  REF == "g"]
dt_sm$chromosome <- chr
sm_out <- long.to.bed(dt_sm[,c("chromosome", "POS")])

# CG bed file
CG_out <- as.data.table(str_locate_all(toupper(paste0(REF, collapse = "")), "CG")[[1]])
CG_out$chromosome <- chr
CG_out <- CG_out[,c("chromosome", "start", "end")]

# bp bed file
bp_out <- data.table(
  chromosome = chr,
  start = start.POS,
  end = end.POS,
  n_N_mask = sum((N_out$end + 1) - N_out$start),
  n_soft_mask = sum((sm_out$end + 1) - sm_out$start),
  n_CG = sum((CG_out$end + 1) - CG_out$start)
)

# convert REF to upper case
ref_out <- data.table(REF = toupper(dt_ref$REF))

### EXPORT

fwrite(N_out, file = out.file.N, col.names = F, append = T, compress = "gzip") # N
fwrite(sm_out, file = out.file.sm, col.names = F, append = T, compress = "gzip") # soft mask
fwrite(CG_out, file = out.file.CG, col.names = F, append = T, compress = "gzip") # CG
fwrite(bp_out, file = out.file.bp, col.names = F, append = T) # bp_counts
fwrite(ref_out, file = out.file.ref, col.names = F, compress = "gzip") # bp_counts



