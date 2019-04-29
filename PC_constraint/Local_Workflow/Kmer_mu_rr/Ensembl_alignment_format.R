# SCRIPT that merges alignment and reference files to include all refernce bases in seqeunce

rm(list=ls())
graphics.off()

library(data.table)


alakazam <- function(in.ref.path, 
                     in.ref.file, 
                     in.align.path, 
                     in.align.file, 
                     chr, 
                     out.path, 
                     out.file){
  for (i in chr){
    
    # read files
    ref.file.path <- paste0(in.ref.path, in.ref.file, i, ".txt")
    ref <- fread(ref.file.path)
    align.file.path <- paste0(in.align.path, in.align.file, i, ".txt")
    align <- fread(align.file.path)
    colnames(align) <- c("POS", "REF", "SPEC")
    
    # merge
    out <- align[ref, on = c("POS", "REF")]
    
    # format
    out$REF[out$REF == "N"] <- NA
    
    # write files
    out.file.path <- paste0(out.path, out.file, i, ".txt")
    fwrite(x = out, file = out.file.path, sep = "\t")
    print (paste("chr", i, "done!", sep = " "))
  }
}

# human <- 
H_in.ref.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_reference/"
H_in.ref.file <- "H_reference_chr"
H_in.align.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_alignments/"
H_in.align.file <- "H_chimpanzee_allignment_chr"
H_chr <- c(1:22, "X")
H_out.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_alignments/"
H_out.file <- "Human_chimpanzee_alignment_chr"

alakazam(in.ref.path = H_in.ref.path, 
          in.ref.file = H_in.ref.file, 
          in.align.path = H_in.align.path, 
          in.align.file = H_in.align.file, 
          chr = H_chr, 
          out.path = H_out.path, 
          out.file = H_out.file)

M_in.ref.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_reference/"
M_in.ref.file <- "M_reference_chr"
M_in.align.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_alignments/"
M_in.align.file <- "M_caroli_allignment_chr"
M_chr <- c(1:19, "X")
M_out.path <- "/Volumes/Untitled/PC_constraint/Paper/Data/Ensembl_alignments/"
M_out.file <- "Mouse_Muscaroli_alignment_chr"

alakazam(in.ref.path = M_in.ref.path, 
         in.ref.file = M_in.ref.file, 
         in.align.path = M_in.align.path, 
         in.align.file = M_in.align.file, 
         chr = M_chr, 
         out.path = M_out.path, 
         out.file = M_out.file)

##### 
# in.ref.path = H_in.ref.path
# in.ref.file = H_in.ref.file
# in.align.path = H_in.align.path
# in.align.file = H_in.align.file
# chr = H_chr
# i = 1
# out.path = H_out.path
# out.file = H_out.file

