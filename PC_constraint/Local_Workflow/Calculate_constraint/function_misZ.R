### FUNCTION that calculates missenseZ and missenseZ percentile
mis_Z <- function(df, maf){
  
  ### FORMAT 
  
  SNVs_cols <-   c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", 
                   "intron_length", "cds_length", "p_syn", "p_mis", "p_non", "r_mu",
                   paste(c("n_synonymous", "n_missense", "n_nonsense", "n_spliceAD", "n_intron"), maf, sep = "_"))
  df_temp <- df[,SNVs_cols, with = F]
  colnames(df_temp) <- c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", 
                         "intron_length", "cds_length", "p_syn", "p_mis", "p_non", "r_mu",
                         "n_syn", "n_mis", "n_non", "n_spliceAD", "n_intron")
  
  ### CREATE VARIABLES 
  
  # calculate intron realive mutation rate
  df_temp$intron_mu <- df_temp$n_intron/df_temp$intron_length
  df_temp$intron_mu[which(df_temp$intron_mu > 1)] <- 1
  
  # calculate "n_functional" and "n_total"
  df_temp$n_total <- df_temp$n_mis + df_temp$n_syn
  
  ### QC
  
  # remove na's
  rm.ind <- which(is.na(df_temp$n_syn) | is.na(df_temp$n_mis) | 
                    is.na(df_temp$p_syn) | is.na(df_temp$p_mis) | 
                    is.na(df_temp$r_mu) | is.na(df_temp$intron_mu))
  if (length(rm.ind) != 0) { df_temp <- df_temp[-rm.ind,] }
  
  # remove genes n_mis or n_syn < 1
  df_temp <- subset(df_temp, df_temp$n_syn >= 1 | df_temp$n_mis >= 1)
  
  # remove genes n_total > 300
  df_temp <- subset(df_temp, df_temp$n_total <= 300)
  
  # remove genes with length > 15000
  df_temp_mod <- subset(df_temp, df_temp$cds_length <= 15000)
  
  # fit model with all genes and remove genes based on Cook's distance 
  # df_mod_all <- lm(df_temp_mod$n_syn ~ df_temp_mod$p_syn + df_temp_mod$r_mu + df_temp_mod$intron_mu)
  # summary(df_mod_all)
  # plot(df_mod_all)
  # df_cooksd <- cooks.distance(df_mod_all)
  # # df_sample_size <- nrow(df_temp)
  # # influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
  # df_influential <- as.numeric(names(df_cooksd)[(df_cooksd > 0.5)])
  # # influential <- as.numeric(names(cooksd)[(cooksd > 4 * mean(cooksd, na.rm = TRUE))])
  # df_temp_removed <- df_temp[df_influential,]
  # df_temp_mod <- df_temp[-df_influential,]
  # 
  ### Calculate missense Z-score
  
  # fit model for synonymous variants
  mod <- lm(n_syn ~ p_syn + r_mu + intron_mu, data = df_temp_mod)
  mse <- mean(mod$residuals^2) # mean squared error
  
  ### PREDICT
  df_temp_mod$pred_mis <- predict(mod, data.frame(p_syn = df_temp_mod$p_mis,
                                                  r_mu = df_temp_mod$r_mu,
                                                  intron_mu = df_temp_mod$intron_mu))
  
  ### missense Z SCORES
  dif <- df_temp_mod$pred_mis - df_temp_mod$n_mis
  mean_dif <- mean(dif, na.rm = T)
  sd_dif <- sd(dif, na.rm = T)
  df_temp_mod$mis_Z <- (dif - mean_dif)/sd_dif
  # df_temp_mod$mis_Z <- scale(dif, center= TRUE, scale=TRUE)
  
  ### Normalise distribution
  # dist <- df_temp_mod$mis_Z[df_temp_mod$mis_Z < 0 & df_temp_mod$mis_Z > -5]
  # dist <- c(dist, dist*-1)
  # sd_dist <- sd(dist)
  # test <- df_temp_mod$mis_Z/sd_dist
  
  # Calculate the percentile of each funZ
  percentile <- ecdf(df_temp_mod$mis_Z)
  df_temp_mod$mis_Z_percentile <- percentile(df_temp_mod$mis_Z)
  
  # subset output
  # df_temp_mod <- df_temp_mod[,c("external_gene_name", "ensembl_gene_id", 
  #                               "ensembl_transcript_id", 
  #                               "mis_Z", "mis_Z_percentile")]
  
  # rename funZ
  names(df_temp_mod)[names(df_temp_mod) == "mis_Z"] <- paste0("mis_Z_", maf)
  names(df_temp_mod)[names(df_temp_mod) == "mis_Z_percentile"] <- paste0("mis_Z_percentile_", maf)
  
  out <- list(df_temp_mod, mod, mse)
  return(out)
}