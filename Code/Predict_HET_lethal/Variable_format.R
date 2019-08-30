rm(list = ls())
graphics.off()

library(data.table)


### IMPORT RAW

# HET failures
hetF <- fread("~/Dropbox/PhD/Data/IMPC/HET_failures_from_Mary.csv")
# IMPC genes
impc <- fread("~/Dropbox/PhD/Data/IMPC/IMPC_ALL_statistical_results.csv")
# Human pLI ans oe and lof_Z
gnomad <- fread("~/Dropbox/PhD/Data/Constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt")
# Mouse funZ 
funZ <- fread("~/Dropbox/PhD/Data/PC_constraint/Constraint/M_funZ_allMUSMUS_v2.csv")
# Human mouse orthologues
HMorth <- fread("~/Dropbox/PhD/Data/Ensembl/BioMart/Ensembl_v94_HM_orthologues.csv")


### FORMAT

## HET failures and IMPC genes

# subset gene names
hetF <- unique(hetF$`Gene Marker Symbol`)
impc <- unique(impc$marker_symbol)
# remove genes from both groups
duplicated_genes <- impc[which(impc %in% hetF)]
hetF <- hetF[-which(hetF %in% duplicated_genes)]
impc <- impc[-which(impc %in% duplicated_genes)]
# add het lethal coliumn and rbind
impc <- rbind(data.table(M_external_gene_name = hetF,
                         Het_fail = rep(1, length(hetF))),
              data.table(M_external_gene_name = impc,
                   Het_fail = rep(0, length(impc)))
)
# remove temp files
rm(hetF, duplicated_genes)

## Orthologues

# subset orthologues with confidence of 1
HMorth <- subset(HMorth, HMorth$orthology_confidence == 1)
# take conservation as smallest value in "Maa_match_Haa" or "Haa_match_Maa"
HMorth <- transform(HMorth, AAconservation = pmin(Haa_match_Maa, Maa_match_Haa))
# subset columns 
HMorth <- HMorth[,c("M_external_gene_name", "H_external_gene_name",
                    "M_ensembl_transcript_id", "H_ensembl_transcript_id",
                    "orthology_type", "AAconservation")]
 
## Mouse funZ

# subset and rename columns 
funZ <- funZ[,c("external_gene_name", "ensembl_transcript_id", "fun_Z")]
colnames(funZ) <- c("M_external_gene_name", "M_ensembl_transcript_id", "fun_Z")


## Human pLI ans oe and lof_Z

# subset and rename columns 
gnomad <- gnomad[,c("gene", "transcript", "oe_lof_upper", "lof_z", "pLI")]
colnames(gnomad) <- c("H_external_gene_name", "H_ensembl_transcript_id", "oe_lof_upper", "lof_z", "pLI")


## Combine

out_file <- impc[HMorth, on = c("M_external_gene_name")]
out_file <- funZ[out_file, on = c("M_external_gene_name", "M_ensembl_transcript_id")]
out_file <- gnomad[out_file, on = c("H_external_gene_name", "H_ensembl_transcript_id")]


### OUTPUT

fwrite(out_file, "~/Dropbox/PhD/Data/Predict_het_fail.csv")


#######





# tran_no <- unique(tran_no[,c("external_gene_name", "ensembl_transcript_id")])
# tran_no <- as.data.table(table(tran_no$external_gene_name))
# impc <- subset(impc, impc$marker_symbol == "Tfeb" | impc$marker_symbol == "Pigq"  | impc$marker_symbol == "Mthfd2" | impc$marker_symbol == "Ascc3") 

