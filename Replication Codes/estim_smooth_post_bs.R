estim_smooth_post_bs <- function(res_equation,coefficients,oprob=FALSE){
  #DESCRIPTION:
  # Evaluates yield/rating equation at any coefficient (can be different from ML coefficients)
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_smooth_post_bs(pre,main,oprob=FALSE)
  #-------------------------------------------------------------------------------
  #INPUT
  # res_equation    results from yield or rating estimation
  # coefficients    median bootstrapped coefficients (or any other alternative set of coefficients)
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
  
  T <- res_equation$T
  k <- res_equation$k
  smooth_coeffs <- res_equation$smooth_coeffs
  lambda <- res_equation$lambda
  Y <- res_equation$pre$Y
  X <- res_equation$pre$X
  if (oprob){
    LL_out <- LL_oprob(Y = Y,X = X,coeffs = coefficients,smooth_coeffs = smooth_coeffs,lambda = lambda,details = TRUE)
    coeffs_null <- 0*coefficients
    coeffs_null[(length(coeffs_null)-1):length(coeffs_null)] <- c(qnorm(mean(Y==1)),qnorm(1-mean(Y==3)))
    LL_null <- LL_oprob(Y = Y,X = X,coeffs = coeffs_null,smooth_coeffs = smooth_coeffs,lambda = lambda,details = TRUE)
    var_data <- 1
    
    R2 <- 1-LL_out$LL_data/LL_null$LL_data
    R2adj <- R2 + (k-2)/LL_null$LL_data            #adjustment: constant is part of the model
    
    AIC <- -2 * LL_out$LL_data + 2*k
    AICc_Breitung <- NULL
    AICc <- -2 * LL_out$LL_data + 2*k + 2*k*(k+1)/(T-k-1)
    BIC <- -2 * LL_out$LL_data + log(T)*k
  }
  else{
    residuals <- Y - X%*%coefficients
    LL_out <- NULL
    LL_out$LL_data <- -T/2 * log(var(residuals[1:T]))
    LL_out$LL_smooth <- sum(log(dnorm(sqrt(lambda) * diff(coefficients[smooth_coeffs]))))
    LL_out$LL <- -LL_out$LL_data - LL_out$LL_smooth
    var_data <- var(residuals[1:T])
    
    R2 <- 1-var_data/var(res_equation$pre$Y)
    R2adj <- R2-(1-R2)*(k-1)/(T-(k-1)-1)            #adjustment: constant is part of the model
    
    AIC <- T*log(var_data) + 2*k
    AICc_Breitung <- log(sum(residuals[1:T]^2))+((2*k+1)/(T-k-2))
    AICc <- T*log(var_data) + 2*k + 2*k*(k+1)/(T-k-1)
    BIC <- T*log(var_data) + k*log(T)
  }
  
  post_bs <- NULL
  post_bs$LL_out <- LL_out
  post_bs$R2 <- R2
  post_bs$R2adj <- R2adj
  post_bs$AIC <- AIC
  post_bs$AICc_Breitung <- AICc_Breitung
  post_bs$AICc <- AICc
  post_bs$BIC <- BIC
  post_bs$var_data <- var_data
  return(post_bs)
}