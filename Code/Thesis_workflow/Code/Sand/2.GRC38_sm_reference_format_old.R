### SCRIPT that formats reference raw fasta files by chromosome
# Output: 
# reference file (cols: REF (all upper))
# soft mask pos (cols: CHR; START; END))
# N mask pos (cols: CHR; START; END))
# bp counts (cols: chromosome; start; end; n_A; n_T; n_C; n_G; n_N_masked; n_soft_masked)

rm(list = ls())
graphics.off()

library(data.table)

### ARGUMENTS

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are arguments: if not, return an error
if (length(args)==0) { stop("Arguments must be supplied", call.=FALSE) } 

# set args variables
in.file <- args[1]
ref.out.file <- args[2]
sm.out.file <- args[3]
N.out.file <- args[4]
bp.out.file <- args[5]
chr <- args[6]
functions.file <- args[7]

# in.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Raw/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.gz"
# ref.out.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_seq_chr1.txt.gz"
# sm.out.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_sm_chr1.bed"
# N.out.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_N_chr1.bed"
# bp.out.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/h.sap_GRCh38_dna_summary.bed"
# chr <- 1
# functions.file <- "/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"

# in.file <- "~/Dropbox/PhD/Data/Ensembl/Reference/Raw/Mus_musculus.GRCm38.dna_sm.chromosome.19.fa"
# ref.out.file <- ""
# sm.out.file <- ""
# N.out.file <- ""
# bp.out.file  <- ""
# chr <- 19
# functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"

### FUNCTIONS
source(functions.file)

### IMPORT
fasta <- fread(paste0("gunzip -cq ", in.file), sep = "\t")
  
### FORMAT
meta <- names(fasta) # get header description
colnames(fasta) <- "V1"

# identify satrt and end POS from header line
meta.split <- strsplit(meta, ":")
start.POS <- as.integer(meta.split[[1]][5])
end.POS <- as.integer(meta.split[[1]][6])

# create vectors for POS and REF
POS <- seq(from = start.POS, to = end.POS)
REF <- paste(fasta$V1, collapse = '')
REF <- strsplit(REF, "")[[1]]

dt_ref <- data.table(POS = POS, 
                     REF = REF)

# check sequence is continuous
is.sequential <- function(x){ all(abs(diff(x)) == 1) }  
if (dt_ref$POS[1]!=1) { stop("Start not equal to one", call.=FALSE) } 
if (!is.sequential(dt_ref$POS)) { stop("POS not sequential", call.=FALSE) } 

# N coordinates bed file
dt_N <- dt_ref[REF == "N" | 
               REF == "n"] # subset N coordinates
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

# convert REF to upper case
ref_out <- data.table(REF = toupper(dt_ref$REF))
  
# bp out file
bp_out <- data.table(
  chromosome = chr,
  start = 1,
  end = nrow(ref_out),
  n_A = length(which(ref_out$REF == "A")),
  n_T = length(which(ref_out$REF == "T")),
  n_C = length(which(ref_out$REF == "C")),
  n_G = length(which(ref_out$REF == "G")),
  n_N_mask = sum((N_out$end + 1) - N_out$start),
  n_soft_mask = sum((sm_out$end + 1) - sm_out$start)
)
  
### EXPORT
  
# ref
fwrite(ref_out, file = ref.out.file, col.names = F, compress = "gzip")
# sm
fwrite(sm_out, file = sm.out.file, col.names = F, compress = "gzip", append = T)
# N
fwrite(N_out, file = N.out.file, col.names = F, compress = "gzip", append = T)
# bp_counts
fwrite(bp_out, file = bp.out.file, col.names = F, append = T)





