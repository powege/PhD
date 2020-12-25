### SCRIPT that formats IMPC raw data.

## Outputs file with columns:
# gene
# viability_status (L, SV; VP; VN)
# fertility_status (FFMF; FIMI; FIMF; FFMI)
# n_parameter_tests
# n_parameter_sig
# parameter_hit_rate
# n_procedure_tests
# n_procedure_sig
# procedure_hit_rate
# n_top_level_MP_tests
# n_top_level_MP
# top_level_MP_hit_rate

## Outputs file with columns:
# gene
# top_level_MP


rm(list = ls())
graphics.off()

library(data.table)
library(tidyr)


### ARGUMENTS

impc.sig.in <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/IMPC/release_12/genotype-phenotype-assertions-IMPC.csv.gz" # significant MP terms
impc.all.in <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/IMPC/release_12/statistical-results-ALL.csv.gz" # all statistical results
out.file.1 <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/IMPC/Formatted/IMPC_r12_viability_fertility_pleiotropy.csv"
out.file.2 <- "~/Dropbox/PhD/Data/Thesis_workflow/Data/IMPC/Formatted/IMPC_r12_top_level_MPs.csv"

### IMPORT

impc.sig <- fread(paste0("gunzip -cq ", impc.sig.in))
impc.all <- fread(paste0("gunzip -cq ", impc.all.in))

### QC

# gene name
impc.sig <- impc.sig[marker_symbol != "" & !is.na(marker_symbol)]
impc.all <- impc.all[marker_symbol != "" & !is.na(marker_symbol)]

# homozygotes
impc.sig <- impc.sig[zygosity == "homozygote"]
impc.all <- impc.all[zygosity == "homozygote"]

# IMPC
impc.sig <- impc.sig[resource_name == "IMPC"]
impc.all <- impc.all[resource_name == "IMPC"]

# filter genes with no significant hits and less than 100 parameter tests
VN_parameter <- unique(impc.all[!impc.all$marker_symbol %in% impc.sig$marker_symbol][,c("marker_symbol", "parameter_stable_id")])
VN_parameter <- as.data.table(table(VN_parameter$marker_symbol))
# hist(VN_parameter$N)
# summary(VN_parameter)
VN_insuf_parameter <- VN_parameter$V1[VN_parameter$N < 100]
impc.all <- impc.all[!marker_symbol %in% VN_insuf_parameter]

### FORMAT 

## FERTILITY

all_tested_male <- impc.all[parameter_name %in% c("Gross Findings Male")] # subset by parameter for testing fertility
all_tested_male <- unique(all_tested_male$marker_symbol)
all_tested_female <- impc.all[parameter_name %in% c("Gross Findings Female")] # subset by parameter for testing fertility
all_tested_female <- unique(all_tested_female$marker_symbol)
all_tested_both <- Reduce(intersect, list(all_tested_male, all_tested_female))
FI <- impc.sig[mp_term_name %in% "female infertility"] # subset by fertility phenotype
FI <- unique(FI$marker_symbol)
MI <- impc.sig[mp_term_name %in% "male infertility"] # subset by fertility phenotype
MI <- unique(MI$marker_symbol)
FIMI <- Reduce(intersect, list(FI, MI))
FI <- FI[!FI %in% FIMI]
MI <- MI[!MI %in% FIMI]
FFMF <- all_tested_both[!all_tested_both %in% c(FIMI, FI, MI)]
Reduce(intersect, list(FIMI, FI, MI, FFMF)) # check groups are mutually exclusive
fertility <- data.table(external_gene_name = c(FIMI, FI, MI, FFMF),
                        fertility_status = c(rep("FIMI", length(FIMI)),
                                             rep("FI", length(FI)),
                                             rep("MI", length(MI)),
                                             rep("FFMF", length(FFMF))))

## VIABILITY

all_ko <- unique(impc.all$marker_symbol)
all_sig <- unique(impc.sig$marker_symbol)
SV <- all_ko[all_ko %in% impc.sig$marker_symbol[impc.sig$mp_term_name == "preweaning lethality, incomplete penetrance"]]# defined by MP of preweaning lethality, incomplete penetrance
L <- all_ko[all_ko %in% impc.sig$marker_symbol[impc.sig$mp_term_name %in% "preweaning lethality, complete penetrance"]] # defined by MP of preweaning lethality, complete penetrance
L <- L[!L %in% SV] # if subviable remove from lethal category
VP <- all_sig[!all_sig %in% c(L, SV)]
VN <- all_ko[!all_ko %in% c(all_sig)]
Reduce(intersect, list(L, SV, VP, VN)) # check groups are mutually exclusive
viability <- data.table(external_gene_name = c(L, SV, VP, VN),
                        viability_status = c(rep("L", length(L)),
                                             rep("SV", length(SV)),
                                             rep("VP", length(VP)),
                                             rep("VN", length(VN))))

## PLEIOTROPY

# Parameters
parameter_all <- unique(impc.all[parameter_stable_id != "" & !is.na(parameter_stable_id)][,c("marker_symbol", "parameter_stable_id")])
parameter_all <- as.data.table(table(parameter_all$marker_symbol))
colnames(parameter_all) <- c("external_gene_name", "n_tests_parameter")
parameter_sig <- unique(impc.sig[parameter_stable_id != "" & !is.na(parameter_stable_id)][,c("marker_symbol", "parameter_stable_id")])
parameter_sig <- as.data.table(table(parameter_sig$marker_symbol))
colnames(parameter_sig) <- c("external_gene_name", "n_sig_parameter")

# Procedure
procedure_all <- unique(impc.all[procedure_stable_id != "" & !is.na(procedure_stable_id)][,c("marker_symbol", "procedure_stable_id")])
procedure_all <- as.data.table(table(procedure_all$marker_symbol))
colnames(procedure_all) <- c("external_gene_name", "n_tests_procedure")
procedure_sig <- unique(impc.sig[procedure_stable_id != "" & !is.na(procedure_stable_id)][,c("marker_symbol", "procedure_stable_id")])
procedure_sig <- as.data.table(table(procedure_sig$marker_symbol))
colnames(procedure_sig) <- c("external_gene_name", "n_sig_procedure")

# Top-level MP
top_level_mp_all <- impc.all[,c("marker_symbol", "top_level_mp_term_name")]
top_level_mp_all <- top_level_mp_all[!is.na(top_level_mp_all$top_level_mp_term_name),]
top_level_mp_all <- top_level_mp_all[top_level_mp_all$top_level_mp_term_name != "",]
top_level_mp_all <- separate_rows(top_level_mp_all, top_level_mp_term_name, sep = ",")
top_level_mp_all <- unique(top_level_mp_all)
top_level_mp_all <- as.data.table(table(top_level_mp_all$marker_symbol))
colnames(top_level_mp_all) <- c("external_gene_name", "n_tests_top_level_mp")
top_level_mp_sig <- impc.sig[, c("marker_symbol", "top_level_mp_term_name")]
colnames(top_level_mp_sig) <- c("external_gene_name", "top_level_mp_term_name")
top_level_mp_sig <- top_level_mp_sig[!is.na(top_level_mp_sig$top_level_mp_term_name),]
top_level_mp_sig <- top_level_mp_sig[top_level_mp_sig$top_level_mp_term_name != "",]
top_level_mp_sig <- separate_rows(top_level_mp_sig, top_level_mp_term_name, sep = ",")
top_level_mp_sig <- unique(top_level_mp_sig)
top_level_mp_sig <- as.data.table(table(top_level_mp_sig$external_gene_name))
colnames(top_level_mp_sig) <- c("external_gene_name", "n_sig_top_level_mp")

# merge
pleiotropy <- parameter_sig[parameter_all, on = "external_gene_name"]
pleiotropy <- procedure_all[pleiotropy, on = "external_gene_name"]
pleiotropy <- procedure_sig[pleiotropy, on = "external_gene_name"]
pleiotropy <- top_level_mp_all[pleiotropy, on = "external_gene_name"]
pleiotropy <- top_level_mp_sig[pleiotropy, on = "external_gene_name"]
pleiotropy[is.na(pleiotropy)] <- 0

# hit rates
pleiotropy$hit_rate_top_level_mp <- pleiotropy$n_sig_top_level_mp / pleiotropy$n_tests_top_level_mp
pleiotropy$hit_rate_procedure <- pleiotropy$n_sig_procedure / pleiotropy$n_tests_procedure
pleiotropy$hit_rate_parameter <- pleiotropy$n_sig_parameter / pleiotropy$n_tests_parameter

# hist(pleiotropy$hit_rate_top_level_mp)
# hist(pleiotropy$hit_rate_procedure)
# hist(pleiotropy$hit_rate_parameter)

## TOP-LEVEL MPs

TLMP<- impc.sig[, c("marker_symbol", "top_level_mp_term_name")]
colnames(TLMP) <- c("external_gene_name", "top_level_mp_term_name")
TLMP <- TLMP[!is.na(TLMP$top_level_mp_term_name),]
TLMP <- TLMP[TLMP$top_level_mp_term_name != "",]
TLMP <- separate_rows(TLMP, top_level_mp_term_name, sep = ",")
TLMP <- unique(TLMP)

### EXPORT

# merge output 1 
out1 <- merge(fertility, viability, all = T)
out1 <- pleiotropy[out1, on = "external_gene_name"]
fwrite(out1, out.file.1)

out2 <- TLMP
fwrite(out2, out.file.2)




###############################################################################################

# x <- impc.all[impc.all$mp_term_name %like% "infert"]
# x <- unique(x[,c("marker_symbol", "mp_term_name")])
# length(unique(x$marker_symbol))

# viability.in <- "~/Dropbox/PhD/Data/IMPC/release_12/viability.csv.gz"
# viability <- fread(paste0("gunzip -cq ", viability.in), fill = F)

# viability <- viability[,c("Gene Symbol", "Viability Phenotype HOMs/HEMIs")]
# colnames(viability) <- c("external_gene_name", "viability")


#####

# # filter genes with no significant parameter tests, that have less than the median number or paramerer tests or procedures in VP group
# # distribution of n_procedure_tests for VP KOs
# VP_procedure <- unique(impc.all[impc.all$marker_symbol %in% VP][,c("marker_symbol", "procedure_stable_id")])
# VP_procedure <- as.data.table(table(VP_procedure$marker_symbol))
# # hist(VP_procedure$N)
# # summary(VP_procedure)
# VP_procedure_median <- median(VP_procedure$N)
# 
# # filter genes with no significant hits and less than 8 procedures (median of VP)
# VN_procedure <- unique(impc.all[!impc.all$marker_symbol %in% impc.sig$marker_symbol][,c("marker_symbol", "procedure_stable_id")])
# VN_procedure <- as.data.table(table(VN_procedure$marker_symbol))
# # hist(VN_procedure$N)
# # summary(VN_procedure)
# VN_insuf_procedure <- VN_procedure$V1[VN_procedure$N < VP_procedure_median]
# 
# # distribution of n_parameter_tests for VP KOs
# VP_parameter <- unique(impc.all[impc.all$marker_symbol %in% VP][,c("marker_symbol", "parameter_stable_id")])
# VP_parameter <- as.data.table(table(VP_parameter$marker_symbol))
# # hist(VP_parameter$N)
# # summary(VP_parameter)
# VP_parameter_median <- median(VP_parameter$N)
# 
# # filter genes with no significant hits and less than 154 parameter tests (median of VP)
# VN_parameter <- unique(impc.all[!impc.all$marker_symbol %in% impc.sig$marker_symbol][,c("marker_symbol", "parameter_stable_id")])
# VN_parameter <- as.data.table(table(VN_parameter$marker_symbol))
# # hist(VN_parameter$N)
# # summary(VN_parameter)
# VN_insuf_parameter <- VN_parameter$V1[VN_parameter$N < VP_parameter_median]

# # identify number of MP tests per KO 
# n.mp.tests <- impc.all[,c("marker_symbol", "mp_term_name")]
# n.mp.tests <- n.mp.tests[!is.na(n.mp.tests$mp_term_name),]
# n.mp.tests <- n.mp.tests[n.mp.tests$mp_term_name != "",]
# n.mp.tests <- separate_rows(n.mp.tests, mp_term_name, sep = ",")
# n.mp.tests <- unique(n.mp.tests)
# n.mp.tests <- as.data.table(table(n.mp.tests$marker_symbol))
# colnames(n.mp.tests) <- c("external_gene_name", "n_MP_tests")
# # n.mp.tests <- n.mp.tests[n_MP_tests >= 10]
# 
# # identify number of top-level MP tests per KO 
# n.tlmp.tests <- impc.all[,c("marker_symbol", "top_level_mp_term_name")]
# n.tlmp.tests <- n.tlmp.tests[!is.na(n.tlmp.tests$top_level_mp_term_name),]
# n.tlmp.tests <- n.tlmp.tests[n.tlmp.tests$top_level_mp_term_name != "",]
# n.tlmp.tests <- separate_rows(n.tlmp.tests, top_level_mp_term_name, sep = ",")
# n.tlmp.tests <- unique(n.tlmp.tests)
# n.tlmp.tests <- as.data.table(table(n.tlmp.tests$marker_symbol))
# colnames(n.tlmp.tests) <- c("external_gene_name", "n_top_level_MP_tests")
# 
# # identify number of unique MPs per KO 
# impc.mp <- impc.sig[, c("marker_symbol", "mp_term_name")]
# colnames(impc.mp) <- c("external_gene_name", "mp_term_name")
# impc.mp <- impc.mp[!is.na(impc.mp$mp_term_name),] # remove genes with no top level phenotype term
# impc.mp <- impc.mp[impc.mp$mp_term_name != "",]
# impc.mp <- separate_rows(impc.mp, mp_term_name, sep = ",") # split MPs into multiple categories
# impc.mp <- unique(impc.mp) # remove duplicates by gene name and mp term
# all.counts <- as.data.frame(table(impc.mp$external_gene_name)) # count terms per gene
# colnames(all.counts) <- c("external_gene_name", "n_MP")
# 
# # identify number of unique top-level MPs per KO 
# impc.tl <- impc.sig[, c("marker_symbol", "top_level_mp_term_name")]
# colnames(impc.tl) <- c("external_gene_name", "top_level_mp_term_name")
# impc.tl <- impc.tl[!is.na(impc.tl$top_level_mp_term_name),]
# impc.tl <- impc.tl[impc.tl$top_level_mp_term_name != "",]
# impc.tl <- separate_rows(impc.tl, top_level_mp_term_name, sep = ",")
# impc.tl <- unique(impc.tl)
# top.counts <- as.data.frame(table(impc.tl$external_gene_name))
# colnames(top.counts) <- c("external_gene_name", "n_top_level_MP")
# 
# ### idenmtify genes with no significant MP
# sig_genes <- unique(impc.sig$marker_symbol)
# all_genes <- unique(impc.all$marker_symbol)
# no_pheno_genes <- all_genes[-which(all_genes %in% sig_genes)]
# no_pheno_genes <- data.frame(external_gene_name = no_pheno_genes,
#                              n_MP = rep(0, length(no_pheno_genes)),
#                              n_top_level_MP = rep(0, length(no_pheno_genes)))

# fertility.in <- "~/Dropbox/PhD/Data/IMPC/release_12/fertility.csv.gz"
# fertility <- fread(paste0("gunzip -cq ", fertility.in), fill = F)
# fertility <- fertility[Zygosity == "homozygote"]

# colnames(fertility)[names(fertility) == "Gene Symbol"] <- "external_gene_name"
# fertility$fertility_status <- NA
# fertility$fertility_status[which(fertility$Sex == "male" & fertility$Phenotype == "Fertile" )] <- "MF"
# fertility$fertility_status[which(fertility$Sex == "male" & fertility$Phenotype == "Infertile" )] <- "MI"
# fertility$fertility_status[which(fertility$Sex == "female" & fertility$Phenotype == "Fertile" )] <- "FF"
# fertility$fertility_status[which(fertility$Sex == "female" & fertility$Phenotype == "Infertile" )] <- "FI"
# m_tmp <- fertility[Sex == "male"][,c("external_gene_name", "fertility_status")]
# f_tmp <- fertility[Sex == "female"][,c("external_gene_name", "fertility_status")]
# m_tmp <- m_tmp[external_gene_name %in% f_tmp$external_gene_name]
# f_tmp <- f_tmp[external_gene_name %in% m_tmp$external_gene_name]
# fertility <- m_tmp[f_tmp, on = "external_gene_name"]
# fertility$fertility_status <- paste0(fertility$fertility_status, fertility$i.fertility_status)
# fertility <- unique(fertility[,c("external_gene_name", "fertility_status")])
# rm(m_tmp, f_tmp)
# # table(fertility$fertility_status)

# x <- impc.sig[, c("marker_symbol", "top_level_mp_term_name", "parameter_stable_id")]
# y <- impc.all[, c("marker_symbol", "top_level_mp_term_name", "parameter_stable_id")]
# x1 <- x[parameter_stable_id == "IMPC_HEM_029_001"]
# y1 <- y[parameter_stable_id == "IMPC_HEM_029_001"]


