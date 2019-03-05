# SCRIPT that compares and plots the relationship between the human and mouse 
# trinucleotide specific probabilites of mutation. 

rm(list = ls())
graphics.off()

library(ggplot2)
library(ggrepel)
library(devtools)

### FUNCTIONS 

### FUNCTION for plotting Pearson's correlation on ggplot
corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}

### FUNCTION for plotting lm equation on ggplot
lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}


### INPUT 

h3.mu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_3mer_mu_rate.table")
m3.mu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_3mer_mu_rate.table")

h5.mu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/H_5mer_mu_rate.table")
m5.mu <- fread("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Data/Mu_rates/M_5mer_mu_rate.table")


### FORMAT 

h3.mu <- h3.mu[,c("k3_from", "k3_to", "k3_mu_rr")]
m3.mu <- m3.mu[,c("k3_from", "k3_to", "k3_mu_rr")]
h5.mu <- h5.mu[,c("k5_from", "k5_to", "k5_mu_rr")]
m5.mu <- m5.mu[,c("k5_from", "k5_to", "k5_mu_rr")]

colnames(h3.mu) <- c("k3_from", "k3_to", "H_k3_mu_rr")
colnames(m3.mu) <- c("k3_from", "k3_to", "M_k3_mu_rr")
colnames(h5.mu) <- c("k5_from", "k5_to", "H_k5_mu_rr")
colnames(m5.mu) <- c("k5_from", "k5_to", "M_k5_mu_rr")

df3 <- h3.mu[m3.mu, on=c("k3_from", "k3_to")]
df5 <- h5.mu[m5.mu, on=c("k5_from", "k5_to")]


### MODEL

mod3 <- lm(df3$H_k3_mu_rr ~ df3$M_k3_mu_rr)
summary(mod3)
cor.test(df3$H_k3_mu_rr, df3$M_k3_mu_rr)

mod5 <- lm(df5$H_k5_mu_rr ~ df5$M_k5_mu_rr)
summary(mod5)
cor.test(df5$H_k5_mu_rr, df5$M_k5_mu_rr)

### PLOT

Name <- paste(df3$k3_from, "to", df3$k3_to, sep = " ")
pearsons3 <- corr_eqn(df3$H_k3_mu_rr, df3$M_k3_mu_rr, digits = 3)
  
p3 <- ggplot(df3, aes(x=M_k3_mu_rr, y=H_k3_mu_rr, label=Name)) +
  geom_point(color = "blue", size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T, color = "red") +
  geom_label_repel(aes(label = ifelse(M_k3_mu_rr>0.1,as.character(Name),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  annotate("text", x = 0.03, y = 0.12, label = pearsons3, colour="black", size = 5, parse=TRUE) +
  # geom_text(aes(x = 0.03, y = 0.12, label = lm_eqn(mod)), parse = TRUE) +
  xlab("Relative mutation rate in the mouse linneage") +
  ylab("Relative mutation rate in the human linneage") +
  theme(axis.text=element_blank(),
  axis.ticks=element_blank()) 
p3

pearsons5 <- corr_eqn(df5$H_k5_mu_rr, df5$M_k5_mu_rr, digits = 3)
p5 <- ggplot(df5, aes(x=M_k5_mu_rr, y=H_k5_mu_rr)) +
  geom_point(color = "blue", size = 2, shape = 1) +
  geom_smooth(method = "lm", se = T, color = "red") +
  # geom_label_repel(aes(label = ifelse(M_k3_mu_rr>0.1,as.character(Name),'')),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50') +
  annotate("text", x = 0.03, y = 0.12, label = pearsons5, colour="black", size = 5, parse=TRUE) +
  # geom_text(aes(x = 0.03, y = 0.12, label = lm_eqn(mod)), parse = TRUE) +
  xlab("Relative mutation rate in the mouse linneage") +
  ylab("Relative mutation rate in the human linneage") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) 
p5


### OUTPUT

ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Supplementary/Figure_k3_cor.jpg", plot = p3, device = "jpeg")
ggsave("~/Dropbox/BitBucket_repos/phd/PC_constraint/Paper/Results/Supplementary/Figure_k5_cor.jpg", plot = p5, device = "jpeg")
