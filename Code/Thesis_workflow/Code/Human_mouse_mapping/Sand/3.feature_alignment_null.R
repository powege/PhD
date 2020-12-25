### TO DO 
# run on rescomp parralell by chromosome
# run for annotations of interest (ie promoter, enhancer, CTCF, TAD, Misc)
# what is the relationship between distance to gene and conservation? Does this need to be accounnted for?

rm(list=ls())
graphics.off()

library(data.table)

# INPUT:
# human annotation bed.gz 
# mouse annotation bed.gz
# human mouse alignment bed.gz
# N masked pos

# OUTPUT:
# colnames: annotation; ind; n_total; n_aligned; n_conserved
  
### SET VARS 

# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there are arguments
# if (length(args)==0) {
#   stop("Arguments must be supplied", call.=FALSE)
# } 
# 
# # set args variables
# h.ann.file <- args[1]
# m.ann.file <- args[2]
# align.file <- args[3]
# N.mask.file <- args[4]
# functions.file <- args[5]
# out.file <- args[6]
# null.size <- args[7]
# chr <- args[8]

h.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/h.sap_GRC38_v101_whole_genome_features_multicell.csv.gz"
m.ann.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/m.mus_GRC38_v101_whole_genome_features_multicell.csv.gz"
align.file <- "~/Dropbox/PhD/Data/Interspecific_SNV_mapping/workflow_v1/formatted/hsap_grch38_v_mmus_grcm38_v101_alignment_chr15.bed.gz"
N.mask.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Ensembl/Reference/Formatted/m.mus_GRCm38_dna_N.bed.gz"
functions.file <- "~/Dropbox/GitHub_repos/Phd/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
out.file <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/Sand/hsap_grch38_v_mmus_grcm38_v101_multicell_feature_alignment_null_local.csv"
null.size <- 100
chr <- 15

# h.ann.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/h.sap_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
# m.ann.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Annotation/Formatted/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell.csv.gz"
# align.file <- "/well/lindgren/George/Data/Thesis_workflow/Data/Ensembl/Alignment/Formatted/hsap_grch38_v_mmus_grcm38_v101_alignment.bed.gz"
# functions.file <- "/well/lindgren/George/Code/Thesis_workflow/Code/Functions/FUNCTIONS.R"
# out.file <- ""
# null.size <- 1
  
### FUNCTIONS
source(functions.file)

### FUNCTION that returns the seqence match between bed file and vector
# (ie for each sequence in bed file1, how many integers overlap with vector?)
bed_match <- function(bed, vec){
  
  require(data.table)
  
  bed_dt <- setDT(bed) # set as bed as data.table
  colnames(bed_dt) <- c("start", "end")
  bed_dt[, ind := .I] # add uniqe index to data.table
  setkey(bed_dt) # sets keys // order data by all columns
  
  vec_dt <- as.data.table(vec, key = 'vec') # convert to data.table
  vec_dt[, vec2 := vec] # dublicate column
  
  # Fast overlap join:
  ans1 = foverlaps(vec_dt, bed_dt, by.x = c('vec', 'vec2'), by.y = c('start', 'end'),
                   type = "within", nomatch = 0L)
  counts <- ans1[, .N, keyby = ind] # count by ind
  # merge to inital data
  bed_dt[, n_match := counts[bed_dt, on = .(ind), x.N]]
  setorder(bed_dt, ind) # reorder by ind to get inital order
  bed_dt[, ind := NULL] # deletes ind colum
  bed_dt[is.na(n_match), n_match := 0L] # NAs is 0 count
  
  return(bed_dt$n_match) 
}

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


### IMPORT 
human_ann <- fread(paste0("gunzip -cq ", h.ann.file))
mouse_ann <- fread(paste0("gunzip -cq ", m.ann.file))
align <- fread(paste0("gunzip -cq ", align.file))
Nmask <- fread(paste0("gunzip -cq ", N.mask.file))

### FORMAT

# define annotation categories of interest 
annotation <- unique(human_ann$category)
annotation <- annotation[which(!annotation %in% c("Intron - distal", "Unannotated"))]

# N mask
colnames(Nmask) <- c("chromosome", "strat", "end")

# subset chromosome
human_ann <- human_ann[chromosome == chr]
Nmask <- Nmask[chromosome == chr]

# QC
human_ann <- human_ann[percentage_Nmask < 1]
mouse_ann <- mouse_ann[percentage_Nmask < 1]

# collaps overlap for huamn annotations 
dt_pop <- collapse.overlap(human_ann[,c("chromosome", "start", "end")])

out_list <- list()
# for loop
  for (i in 1:null.size){
    
    # set output vecs
    n_total <- rep(NA, length(annotation))
    n_aligned <- rep(NA, length(annotation))
    n_conserved <- rep(NA, length(annotation))
    
    for (ann.ind in 1:length(annotation)){
    
    human_sub <- human_ann[category == annotation[ann.ind]]
      
    # get a vector of sequence lengths 
    human_sub[, length := 1 + end-start]
    lengths <- human_sub$length
    
    # sample positions from population (all pos)
    sample <- bed.sample.pos(bed = dt_pop, nPOS = length(lengths))
    
    # subtract lengths from end 
    sample$length <- lengths
    sample <- sample[, start := ((start + 1) - length)][,c("chromosome", "start", "end")]
    
    # calculate Nmask and filter Nmask > 1%
    sample$n_Nm <- bed_match(bed = sample[,c("start", "end")],
                                 vec = unlist(seq2(from = Nmask$strat, to = Nmask$end)))
    sample$percentage_Nm <- (sample$n_Nm / (((sample$end +1) - sample$start))) * 100 
    sample <- sample[percentage_Nm < 1]
    
    # calculate total bases
    sample_colapse <- collapse.overlap(sample[, c("chromosome", "start", "end")])
    n_total[ann.ind] <- sum(abs( (sample_colapse$end + 1) - sample_colapse$start))
    
    # orthologous sequences for human annotation to get total alignment
    sample_align <- orthologous.seq(alignment = align[,c("chromosome_human", "start_human", "end_human",
                                                         "chromosome_mouse", "start_mouse", "end_mouse")],
                                    bed = sample_colapse)
    n_aligned[ann.ind] <- sum(abs( (sample_align$end_A + 1) - sample_align$start_A))
    
    # total of mouse annotation in orthologous sequences for conservation
    sample_conserved <- bed.intersect(bed1 = sample_align[,c("chromosome_B", "start_B", "end_B")],
                                      bed2 = mouse_ann[category == "Promoter"][,c("chromosome", "start", "end")])
    n_conserved[ann.ind] <- sum(abs ((sample_conserved$end + 1) - sample_conserved$start))
    
    print(annotation[ann.ind])
  }
  output <- data.table(category = annotation,
                    chromosome = chr,
                                    ind = i,
                                    n_total = n_total,
                                    n_aligned = n_aligned,
                                    n_conserved = n_conserved)
  
  out_list[[i]] <- output
  
  ### EXPORT
  fwrite(output, out.file, append = T)
  
  print(i)
}
out_dt <- do.call("rbind", out_list)
# out_dt$percent_aligned <- ( out_dt$n_aligned / out_dt$n_total ) * 100
# out_dt$percent_conserved <- ( out_dt$n_conserved / out_dt$n_total ) * 100
# hist(out_dt$percent_aligned)
# hist(out_dt$percent_conserved)





