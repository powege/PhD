# SCRIPT that converts fasta format to text file with each base and genomic coordinates
# ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/README

rm(list = ls())
graphics.off()

library(data.table)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Exactly three arguments must be supplied", call.=FALSE)
} 

# set args variables
in.file <- args[1]
out.file <- args[2]
chr <- args[3]

### FUNCTIONS

fasta_convert <- function(in.file, out.file, chr){
    
    # read file
    in.file.path <- paste0(in.file, chr, ".fa")
    fasta <- fread(in.file.path, sep = "\t")
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
    
    out <- data.table(POS = POS, 
                      REF = REF)
    
    # identify masked POS
    ind <- which(         out$REF == "N" |
                          out$REF == "n" | 
                           out$REF == "a" | 
                           out$REF == "c" | 
                           out$REF == "t" |
                           out$REF == "g")
    out$REP_MASK <- 0
    out$REP_MASK[ind] <- 1
    
    # convert REF to upper case
    out$REF <- toupper(out$REF)
    
    out.file.path <- paste0(out.file, chr, ".txt")
    fwrite(out, file = out.file.path, sep = "\t")
    print (paste(chr, "done!", sep = " "))
  
}

### RUN
fasta_convert(in.file, out.file, chr)


