rm(list = ls())
graphics.off()

set.seed(1)

library('dismo')
library('gbm')
library('readr')
library(pROC)
library(plyr)
library(data.table)


##############
### SCRIPT ###
##############

### IMPORT
dt <- fread("~/Dropbox/PhD/Data/Predict_het_fail.csv")

### FORMAT

# subset impc genes for train/test
mod_var <- dt[!is.na(dt$Het_fail),]

# subset one-to-one orthologues
mod_var <- subset(mod_var, mod_var$orthology_type == "ortholog_one2one")
# mod_var$orthology_type <- as.factor(mod_var$orthology_type)

### FIT BRT

model <- gbm.step(data=mod_var,
                  gbm.x = c("fun_Z",
                            "oe_lof_upper",
                            # "lof_z",
                            "pLI",
                            # "orthology_type",
                            "AAconservation"),
                  gbm.y = "Het_fail",
                  n.folds = 10,
                  family = "bernoulli",
                  tree.complexity = 5,
                  learning.rate = 0.005,
                  var.monotone = c(1,1,1,1),
                  bag.fraction = 0.5
)

### EXTRACT MODEL PARAMETERS

# summary
summary(model)

# mean CV ROC
mean_ROC <- mean(model$cv.roc.matrix)

# plot 
gbm.plot(model,
         # n.plots=3, 
         write.title = FALSE,
         plot.layout=c(2, 2)
         )

# variable rel_inf
rel_inf <- relative.influence(model)
rel_inf <- (rel_inf/sum(rel_inf))*100
