rm(list = ls())
graphics.off()

# x <- fread("~/Dropbox/PhD/Data/NC_constraint/Figures_and_tables/Raw/ClinVar_SNVs_by_percentile_650_50.csv")
# table(x$Annotation)
# sum(x$N_ClinVar[x$Annotation == "Human ClinVar all"])
# sum(x$N_ClinVar[x$Annotation == "Human ClinVar aligned"])
# sum(x$N_ClinVar[x$Annotation == "Mouse ClinVar aligned"])
# plot(x$N_ClinVar[x$Annotation == "Human ClinVar all"]~
#        x$Percentile[x$Annotation == "Human ClinVar all"])
# plot(x$N_ClinVar[x$Annotation == "Human ClinVar aligned"]~
#        x$Percentile[x$Annotation == "Human ClinVar aligned"])
# plot(x$N_ClinVar[x$Annotation == "Mouse ClinVar aligned"]~
#        x$Percentile[x$Annotation == "Mouse ClinVar aligned"])


cv <- fread("~/Dropbox/PhD/Data/ClinVar/formatted/ClinVar_human_to_mouse_alignment.txt")
H_con <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_windowhuman_650_50.csv")
M_con <- fread("~/Dropbox/PhD/Data/NC_constraint/Constraint/Constraint_by_windowmouse_650_50.csv")


H_con <- H_con[,c("CHR", "POS_from", "POS_to", "Constraint_percentile_CpG")]
H_con$CHR <- as.integer(H_con$CHR)
M_con <- M_con[,c("CHR", "POS_from", "POS_to", "Constraint_percentile_CpG")]
M_con$CHR <- as.integer(M_con$CHR)

cv <- subset(cv, cv$M_CHR %in% c(1:19))
m_cv <- cv[,c("M_CHR", "M_POS", "M_REF")]
colnames(m_cv) <- c("CHR", "POS", "REF")
m_cv$CHR <- as.integer(m_cv$CHR)
h_cv <- cv[,c("CHR", "POS", "REF")]
h_cv$CHR <- as.integer(h_cv$CHR)

H_out <- setDT(h_cv)[H_con, RANK := Constraint_percentile_CpG, on = .(CHR, POS >= POS_from, POS <= POS_to)]
M_out <- setDT(m_cv)[M_con, RANK := Constraint_percentile_CpG, on = .(CHR, POS >= POS_from, POS <= POS_to)]




df <- data.frame(H = H_out$RANK, M = M_out$RANK)
df <- df[complete.cases(df),]
plot(df$H, df$M)
cor.test(df$H, df$M)
t.test(df$H, df$M, paired = T)
mod <- lm(df$H~df$M)
summary(mod)



#### STACK OVERFLOW

set.seed(1)

df1 <- data.frame(CAT = c(rep(1, 1000000), rep(2, 1000000), rep(3, 1000000)),
                  START = c(seq(1, 10000000, 10), seq(1, 10000000, 10), seq(1, 10000000, 10)),
                  END = c(seq(10, 10000000, 10), seq(10, 10000000, 10), seq(10, 10000000, 10)),
                  RANK = sample(1:100, 3000000, replace = T))

df2 <- data.frame(CAT = sample(1:3, 100, replace = T),
                  INT = sample(1:15000000, 100))

system.time({ 
  
  out <- rep(NA, nrow(df2))
  for (i in 1:nrow(df2)){
    x <- subset(df1, df1$CAT == df2$CAT[i] &
                  df1$START <= df2$INT[i] & 
                  df1$END >= df2$INT[i])
    if (nrow(x) != 0){ out[i] <- x$RANK[1] }
    print(i)
  }
  
})

test <- setDT(df2)[df1, RANK := RANK, on = .(CAT, INT >= START, INT <= END)]
