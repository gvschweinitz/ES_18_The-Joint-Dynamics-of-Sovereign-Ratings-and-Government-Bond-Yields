estim_smooth_post <- function(pre,main,oprob=FALSE){
  # DESCRIPTION:
  # Performs post estimation for a linear model or an ordered probit
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_smooth_post(pre,main,oprob=FALSE)
  #-------------------------------------------------------------------------------
  #ARGUMENTS
  # pre             result from pre-estimation (estimate_smooth_pre.R)
  # main            result from main estimation (estimate_smooth_main.R)
  # oprob           boolean to indicate ordered probit
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
  
  res_smooth <- main$res_smooth
  T <- main$T
  k <- main$k
  smooth_coeffs <- main$smooth_coeffs
  coeffs <- main$res_smooth$coefficients
  lambda <- main$lambda
  
  if (oprob){
    LL_out <- LL_oprob(Y = pre$Y,X = pre$X,coeffs = coeffs,smooth_coeffs = smooth_coeffs,lambda = lambda,details = TRUE)
    var_data <- 1
    
    AIC <- -2 * LL_out$LL_data + 2*k
    AICc_Breitung <- NULL
    AICc <- -2 * LL_out$LL_data + 2*k + 2*k*(k+1)/(T-k-1)
    BIC <- -2 * LL_out$LL_data + log(T)*k
  }
  else{
    LL_out <- NULL
    LL_out$LL_data <- -T/2 * log(var(res_smooth$resid[1:T]))
    LL_out$LL_smooth <- sum(log(dnorm(sqrt(lambda) * diff(res_smooth$coefficients[smooth_coeffs]))))
    LL_out$LL <- -LL_out$LL_data - LL_out$LL_smooth
    var_data <- var(res_smooth$residuals[1:T])
    
    AIC <- T*log(var_data) + 2*k
    AICc_Breitung <- log(sum(res_smooth$residuals[1:T]^2))+((2*k+1)/(T-k-2))
    AICc <- T*log(var_data) + 2*k + 2*k*(k+1)/(T-k-1)
    BIC <- T*log(var_data) + k*log(T)
  }
  
  post <- NULL
  post$LL_out <- LL_out
  post$AIC <- AIC
  post$AICc_Breitung <- AICc_Breitung
  post$AICc <- AICc
  post$BIC <- BIC
  post$var_data <- var_data
  return(post)
}