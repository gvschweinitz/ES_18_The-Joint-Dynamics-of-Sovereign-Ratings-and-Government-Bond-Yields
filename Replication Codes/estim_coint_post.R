estim_coint_post <- function(pre,main){
  # DESCRIPTION:
  # Performs post estimation for a cointegration model
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_coint_post(pre,main)
  #-------------------------------------------------------------------------------
  #INPUT
  # pre             result from pre-estimation (estimate_coint_pre.R)
  # main            result from main estimation (estimate_coint_main.R)
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   post$LL_out   list with likelihood entries. Namely NEGATIVE LL, LL_data and LL_smooth
  #   post$AIC      Akaike
  #   post$AICc_Breitung  Akaike in the way of Breitung, if available
  #   post$AICc     corrected Akaike
  #   post$BIC      BIC
  #   post$var_data variance of residuals
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  T <- 2*main$T
  k <- main$k
  LL_data <- main$detailed_res$LL_lin + main$detailed_res$LL_oprob
  
  LL_out <- NULL
  LL_out$LL <- main$detailed_res$LL
  LL_out$LL_lin <- main$detailed_res$LL_lin
  LL_out$LL_oprob <- main$detailed_res$LL_oprob
  LL_out$LL_data <- LL_data
  LL_out$LL_smooth_lin <- main$detailed_res$LL_smooth_lin
  LL_out$LL_smooth_oprob <- main$detailed_res$LL_smooth_oprob
  LL_out$LL_smooth <- LL_out$LL_smooth_lin + LL_out$LL_smooth_oprob
  
  AIC <- -2 * LL_data + 2*k
  AICc_Breitung <- NULL
  AICc <- -2 * LL_data + 2*k + 2*k*(k+1)/(T-k-1)
  BIC <- -2 * LL_data + log(2*T)*k
  
  post <- NULL
  post$LL_out <- LL_out
  post$AIC <- AIC
  post$AICc_Breitung <- AICc_Breitung
  post$AICc <- AICc
  post$BIC <- BIC
  return(post)
}