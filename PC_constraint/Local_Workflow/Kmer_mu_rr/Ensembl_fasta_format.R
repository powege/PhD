# SCRIPT that converts fasta format to text file with each base and genomic coordinates
# ftp://ftp.ensembl.org/pub/release-94/fasta/mus_musculus/dna/README

rm(list = ls())
graphics.off()

library(data.table)

fasta_convert <- function(in.path, in.file, chr, out.path, out.file){
  for (i in chr){
    
    # read file
    in.file.path <- paste0(in.path, in.file, i, ".fa")
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
    
    out.file.path <- paste0(out.path, out.file, i, ".txt")
    fwrite(out, file = out.file.path, sep = "\t")
    print (paste("chr", i, "done!", sep = " "))
  }
}

# Human
H_in.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_reference/"
H_in.file <- "Homo_sapiens.GRCh38.dna.chromosome."
H_chr <- c(1:22, "X")
H_out.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_reference/"
H_out.file <- "H_reference_chr"

fasta_convert(H_in.path, H_in.file, H_chr, H_out.path, H_out.file)

# Mouse
M_in.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_reference/"
M_in.file <- "Mus_musculus.GRCm38.dna.chromosome."
M_chr <- c(1:19, "X")
M_out.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_reference/"
M_out.file <- "M_reference_chr"

fasta_convert(M_in.path, M_in.file, M_chr, M_out.path, M_out.file)


####

# in.path <- H_in.path
# in.file <- H_in.file
# chr <- c(1:3)
# out.path <- H_out.path
# out.file <- H_out.file
# i <- 1


