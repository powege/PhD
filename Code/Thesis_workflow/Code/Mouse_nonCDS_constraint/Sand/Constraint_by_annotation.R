rm(list = ls())
graphics.off()

library(data.table)
library(MASS)
library(plyr)

## QC:
#   coverage mask
#   soft mask
#   N mask

## Constraint method 1:
#
# train model:
#   create training dataset by sampling sequences (same length, X times as many) from intron distal and unannotated regions.
#   or
#   use intron distal and unannotated sequences as training data 
#   build model to predict n snv ~ covariates (try lm and brt)
#     
# test model: 
#   for each annotaiton:
#     use training model to predict expected n snv ~ covariates
#   
# calculate constraint:
#   calculate constraintn for each annotation as O/E 
#   use prediction interval to define confidence around E
#
## Constraint method 2:
# for each annotation:
#   regress n snvs ~ covariates
#   take residual as measure of constraint

### IMPORT 
ann <- fread("gunzip -cq ~/Dropbox/PhD/Data/Thesis_workflow/Data/Mouse_nonCDS_constraint/Sand/m.mus_GRC38_Ensembl_v101_whole_genome_features_multicell_covariates.csv.gz")

### FORMAT

ann$percentage_DPmask <- ( ann$n_DPmask / ann$length ) * 100

### QC 

ann <- ann[percentage_Nmask < 1]
ann <- ann[percentage_sm < 90]
ann <- ann[percentage_DPmask < 50]

### METHOD 1

train <- ann[category %in% c("Intron - distal", "Unannotated")]
hist(train$length)

mod_train <- lm(formula = n_snv ~ length + n_DPmask + n_soft_mask + n_meth_CG,
                data = train)
summary(mod_train)
# plot(mod_train)

ann_sub <- ann[category == "Exon - UTR"]
obs_exp <- function(ann_sub){
  ann_pred <- as.data.table(predict(mod_train, 
                                         ann_sub[,c("n_snv", "length", "n_DPmask", "n_soft_mask", "n_meth_CG")], 
                                         interval="predict", 
                                         level=0.95))
  colnames(ann_pred) <- c("exp_snv", "exp_snv_lwr", "exp_snv_upr")
  ann_sub <- cbind(ann_sub, ann_pred)
  ann_sub$exp_snv_lwr[which(ann_sub$exp_snv_lwr < 0)] <- 0
  ann_sub$exp_snv[which(ann_sub$exp_snv < 0)] <- 0
  ann_sub$exp_snv_upr[which(ann_sub$exp_snv_upr < 0)] <- 0
  ann_sub$OE_snv <- ann_sub$n_snv/ann_sub$exp_snv
  ann_sub$OE_snv_lwr <- ann_sub$n_snv/ann_sub$exp_snv_lwr
  ann_sub$OE_snv_upr <- ann_sub$n_snv/ann_sub$exp_snv_upr
  ann_sub$OE_snv[which(ann_sub$n_snv == 0 & ann_sub$exp_snv == 0)] <- 1
  ann_sub$OE_snv_lwr[which(ann_sub$n_snv == 0 & ann_sub$exp_snv_lwr == 0)] <- 1
  ann_sub$OE_snv_upr[which(ann_sub$n_snv == 0 & ann_sub$exp_snv_upr == 0)] <- 1
  # hist(ann_sub$OE_snv)
  # hist(ann_sub$OE_snv_upr)
  return(ann_sub)
}

test <- ann[category %in% c("Exon - UTR", 
                            "Promoter", 
                            "Enhancer - proximal", 
                            "Enhancer - distal", 
                            "CTCF binding",        
                            "Miscellaneous",      
                            "TAD boundry")]
oe_output <- ddply(test, "category", obs_exp)


# ### METHOD 2
# 
# # hist(promoter$length)
# mod_promoter_all <- lm(formula = n_snv ~ length + n_DPmask + n_soft_mask + n_meth_CG,
#                 data = promoter)
# summary(mod_promoter_all)
# # plot(mod_promoter_all)
# promoter$RVIS_specific <- studres(mod_promoter_all)
# 
# # cooksd <- cooks.distance(mod_promoter_all)
# # sample_size <- nrow(promoter)
# # influential_ind <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
# # mod_promoter_cook <- lm(formula = n_snv ~ length + n_DPmask + n_soft_mask + n_meth_CG,
# #                         data = promoter[-influential_ind,])
# # summary(mod_promoter_all)
# # plot(mod_promoter_all)
# 
# # hist(utr$length)
# mod_utr_all <- lm(formula = n_snv ~ length + n_DPmask + n_soft_mask + n_meth_CG,
#                    data = utr)
# summary(mod_utr_all)
# # plot(mod_utr_all)
# utr$RVIS_specific <- studres(mod_utr_all)
# 
# output <- rbind(promoter, utr, fill = T)

### EXPORT

fwrite(oe_output, "~/Dropbox/PhD/Data/Thesis_workflow/Data/Sand/m.mus_GRC38_Ensembl_v101_regulatory_features_multicell_constraint.csv.gz")


  