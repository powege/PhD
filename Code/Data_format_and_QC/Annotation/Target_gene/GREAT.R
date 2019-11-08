rm(list = ls())
graphics.off()

library(data.table)
library(rGREAT)

### FUNCTIONS

alakazam <- function(input, species, version, rule){
  
  job = submitGreatJob(input, 
                       species = species,  
                       version = version,
                       rule = rule,
                       includeCuratedRegDoms = FALSE)
  
  res = plotRegionGeneAssociationGraphs(job, type = 2)
  
  out <- data.frame(chromosome = seqnames(res),
                    start = start(res),
                    end = end(res),
                    gene = elementMetadata(res)$gene,
                    TSS_dist = elementMetadata(res)$distTSS
  )
  out <- setDT(out)
  
  out <- out[input, on = c("chromosome", "start", "end")]
  
  return(out) 
}

### SET ARGS

bed_file <- "~/Dropbox/PhD/Data/Ensembl/Annotation/Human_GRC38_GENCODE_RegBuild_annotation_seqQC.csv"

### IMPORT

bed <- fread(bed_file)

### FORMAT

colnames(bed) <- c("chromosome", "start", "end", "category", "strand", "length", "f_mask", "f_CG")
bed$chromosome <- paste0("chr", bed$chromosome)
input_promoter <- bed[category == "Promoter"]
input_other <- bed[category %in% c("CTCF binding", "Enhancer - distal", "Enhancer - proximal")]

output_promoter <- alakazam(input = input_promoter, 
                 species = "hg38",
                 version = "4.0.4",
                 rule = "oneClosest")

# out_list <- list()
# for (chr in 1:22){
#   input_sub <- input_promoter[chromosome == paste0("chr", chr)]
#   out_list[[chr]] <- alakazam(input = input_sub, 
#                    species = "hg38", 
#                    version = "4", 
#                    rule = "oneClosest")
# }








