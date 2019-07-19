## FUNCTION that calculates funZ and funZ percentile
fun_Z <- function(df, fun_maf, syn_maf){
  
  ### FORMAT 
  SNVs_cols <-   c("chromosome_name", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", 
                   "intron_length", "cds_length", "p_syn", "p_mis", "p_non", "r_mu",
                   paste(c("n_synonymous"), syn_maf, sep = "_"),
                   paste(c("n_missense", "n_nonsense", "n_spliceAD", "n_intron"), fun_maf, sep = "_"))
  df_temp <- df[,SNVs_cols, with = F]
  colnames(df_temp) <- c("chromosome_name", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", 
                         "intron_length", "cds_length", "p_syn", "p_mis", "p_non", "r_mu",
                         "n_syn", "n_mis", "n_non", "n_spliceAD", "n_intron")
  
  ### CREATE VARIABLES 
  
  # calculate intron realive mutation rate
  df_temp$intron_mu <- df_temp$n_intron/df_temp$intron_length
  df_temp$intron_mu[which(df_temp$intron_mu > 1)] <- 1
  
  # calculate "n_functional" and "n_total"
  df_temp$n_functional <- df_temp$n_mis + df_temp$n_non 
  df_temp$n_total <- df_temp$n_functional + df_temp$n_syn
  
  # calculate "p_functional" and "p_total"
  df_temp$p_functional <- df_temp$p_mis + df_temp$p_non
  df_temp$p_total <- df_temp$p_functional + df_temp$p_syn

  # calculate expected proportion of functional and synonymous
  df_temp$exprop_syn <- df_temp$n_total*(df_temp$p_syn/df_temp$p_total)
  df_temp$exprop_fun <- df_temp$n_total*(df_temp$p_functional/df_temp$p_total)
  
  ### QC
  
  # remove na's
  rm.ind <- which(is.na(df_temp$n_syn) | is.na(df_temp$n_total) |  is.na(df_temp$p_syn) |
                  is.na(df_temp$p_functional) |  is.na(df_temp$r_mu) | is.na(df_temp$intron_mu))
  if (length(rm.ind) != 0) { df_temp <- df_temp[-rm.ind,] }
  
  # remove genes n_total < 1
  df_temp <- subset(df_temp, df_temp$n_total >= 1)
  
  # remove genes n_total > 300
  df_temp <- subset(df_temp, df_temp$n_total <= 300)
  
  # remove genes with length > 15000
  df_temp_mod <- subset(df_temp, df_temp$cds_length <= 15000)
  
  # # fit model with all genes and remove genes based on Cook's distance
  # df_mod_all <- lm(df_temp$n_syn ~ df_temp$n_total + df_temp$p_syn + df_temp$r_mu + df_temp$intron_mu)
  # summary(df_mod_all)
  # plot(df_mod_all)
  # df_cooksd <- cooks.distance(df_mod_all)
  # # df_sample_size <- nrow(df_temp)
  # # influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
  # df_influential <- as.numeric(names(df_cooksd)[(df_cooksd > 0.5)])
  # # influential <- as.numeric(names(cooksd)[(cooksd > 4 * mean(cooksd, na.rm = TRUE))])
  # df_temp_removed <- df_temp[df_influential,]
  # df_temp_mod <- df_temp[-df_influential,]
  
  ### Calculate Functional Z-score
  
  # fit model for synonymous variants
  mod <- lm(n_syn ~ 
              p_syn +
              # n_total +
              # exprop_syn +
              r_mu + 
              intron_mu, data = df_temp_mod)
  # plot(mod)
  # summary(mod)
  
  # model coeff
  coef <- summary(mod)$coef
  r2 <- summary(mod)$adj.r.squared
  mse <- mean(mod$residuals^2) # mean squared error
  
  ### PREDICT
  df_temp_mod$exp_fun <- predict(mod, data.frame(n_total = df_temp_mod$n_total, 
                                                 exprop_syn = df_temp_mod$exprop_fun,
                                                  p_syn = df_temp_mod$p_functional,
                                                  r_mu = df_temp_mod$r_mu,
                                                  intron_mu = df_temp_mod$intron_mu))
  
  ### funZ SCORES
  dif <- df_temp_mod$exp_fun - df_temp_mod$n_functional
  mean_dif <- mean(dif, na.rm = T)
  sd_dif <- sd(dif, na.rm = T)
  df_temp_mod$fun_Z <- (dif - mean_dif)/sd_dif
  # df_temp_mod$fun_Z <- scale(dif, center= TRUE, scale=TRUE)
  
  ### Normalise distribution
  # dist <- df_temp_mod$fun_Z[df_temp_mod$fun_Z < 0 & df_temp_mod$fun_Z > -5]
  # dist <- c(dist, dist*-1)
  # sd_dist <- sd(dist)
  # test <- df_temp_mod$fun_Z/sd_dist
  
  # Calculate the percentile of each funZ
  percentile <- ecdf(df_temp_mod$fun_Z)
  df_temp_mod$fun_Z_percentile <- percentile(df_temp_mod$fun_Z)
  
  # subset output
  # df_temp_mod <- df_temp_mod[,c("chromosome_name", "external_gene_name", "ensembl_gene_id", 
  #                               "ensembl_transcript_id",
  #                               "fun_Z", "fun_Z_percentile")]
  
  # rename funZ
  names(df_temp_mod)[names(df_temp_mod) == "fun_Z"] <- paste0("fun_Z_", fun_maf)
  names(df_temp_mod)[names(df_temp_mod) == "fun_Z_percentile"] <- paste0("fun_Z_percentile_", fun_maf)
  
  out <- list(df_temp_mod, mod, mse)
  return(out)
}