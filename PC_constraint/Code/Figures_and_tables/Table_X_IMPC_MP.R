rm(list = ls())
graphics.off()

source("~/Dropbox/GitHub_repos/PhD/PC_constraint/Code/Figures_and_tables/Figure_X_IMPC_MP.R")
library(gridExtra)

### Wilcoxon
# test FunZ for differneces in funZ between MP groupings
terms <- levels(temp$top_level_mp_term_name)
N_sig <- rep(NA, length(terms))
N_nosig <- rep(NA, length(terms))
impc.statistic.f <- rep(NA, length(terms))
impc.p.val.f <- rep(NA, length(terms))
mean_sig <- rep(NA, length(terms))
mean_nosig <- rep(NA, length(terms))

for (i in 1:length(terms)){
  
  N_sig[i] <- length(temp$gene_symbol[temp$top_level_mp_term_name == terms[i] & temp$sig == 1])
  N_nosig[i] <- length(temp$gene_symbol[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
  
  impc.test <- wilcox.test(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 1], 
                           temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
  impc.statistic.f[i] <- impc.test$statistic
  impc.p.val.f[i] <- impc.test$p.value
  mean_sig[i] <- mean(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 1])
  mean_nosig[i] <- mean(temp$funZ_IMPC_percentile[temp$top_level_mp_term_name == terms[i] & temp$sig == 0])
}
out <- data.frame(
                  top_level_MP_term = terms,
                  n_KO_with_MP = N_sig,
                  n_KO_without_MP = N_nosig,
                  mean_funZ_with_MP = round(mean_sig, digits = 3),
                  mean_funZ_without_MP = round(mean_nosig, digits = 3),
                  statistic = formatC(impc.statistic.f, format = "e", digits = 2),
                  p_val = formatC(p.adjust(impc.p.val.f, method = "bonferroni", n = length(impc.p.val.f)), format = "e", digits = 2)
                  )


fwrite(out, "~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/Table_X_IMPC_MP.txt")
jpeg("~/Dropbox/PhD/Data/PC_constraint/Figures_and_tables/Table_X_IMPC_MP.jpeg", height = 22*nrow(out), width = 125*ncol(out))
grid.table(out, rows = NULL)
dev.off()




