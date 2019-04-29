### FUNCTION that calculates RVIS and RVIS percentile
# human functional MAF > 0.0005, total MAF > 0.0001 (all variants)
RVIS <- function(df, fun_maf, all_maf){
  
  ### FORMAT 
  
  SNVs_cols <-   c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
                   paste(c("n_synonymous", "n_missense", "n_nonsense", "n_spliceAD"), all_maf, sep = "_"),
                   paste(c("n_missense", "n_nonsense", "n_spliceAD"), fun_maf, sep = "_"))
  df_temp <- df[,SNVs_cols, with = F]
  colnames(df_temp) <- c("external_gene_name","ensembl_gene_id","ensembl_transcript_id",        
                         "n_synonymous_all","n_missense_all","n_nonsense_all","n_spliceAD_all",
                         "n_missense_maf", "n_nonsense_maf", "n_spliceAD_maf")
  
  ### CREATE VARIABLES 
  
  # calculate "n_functional" and "n_total"
  df_temp$n_functional <- df_temp$n_missense_maf + df_temp$n_nonsense_maf + df_temp$n_spliceAD_maf
  df_temp$n_total <- df_temp$n_missense_all + df_temp$n_nonsense_all + df_temp$n_spliceAD_all + df_temp$n_synonymous_all
  
  ### QC
  
  # remove na's
  rm.ind <- which(is.na(df_temp$n_total) | is.na(df_temp$n_functional))
  if (length(rm.ind) != 0) { df_temp <- df_temp[-rm.ind,] }
  
  # remove genes n_total < 1
  df_temp <- subset(df_temp, df_temp$n_total >= 1)
  
  ## Fit linear model without weights 
  mod <- lm(n_functional ~ n_total, data = df_temp)
  # plot(n_functional ~ n_total, data = df_temp)
  # summary(mod)
  # plot(mod)
  
  # Take RVIS as studentised residual
  df_temp$RVIS <- studres(mod)
  
  # invert scores 
  df_temp$RVIS <- 0-df_temp$RVIS
  
  # Calculate the parcentile of each RVIS
  percentile <- ecdf(df_temp$RVIS)
  df_temp$RVIS_percentile <- percentile(df_temp$RVIS)
  
  # subset output
  # df_temp <- df_temp[,c("external_gene_name", "ensembl_gene_id", 
  #                               "ensembl_transcript_id", 
  #                       "RVIS", "RVIS_percentile")]
  
  # rename funZ
  names(df_temp)[names(df_temp) == "RVIS"] <- paste0("RVIS_", fun_maf)
  names(df_temp)[names(df_temp) == "RVIS_percentile"] <- paste0("RVIS_percentile_", fun_maf)
  
  out <- list(df_temp, mod)
  
  return(out)
}
