LL_oprob <- function(Y,X,coeffs,smooth_coeffs,lambda,details=FALSE) {
  # DESCRIPTION
  # This function calculates the log-likelihood for a ordered probit with (partially) smoothed coefficients.
  # The log-likelihood has the components:
  #     - fitting part, summed over observations
  #     - smoothness at (differenced) smoothed coefficients
  #----------------------------------------------------
  # USAGE
  # LL_oprob(Y,X,coeffs,smooth_coeffs,lambda,details = FALSE)
  #----------------------------------------------------
  #INPUT
  # Y             vector of the dependent variable (discrete number of values == classes)
  # X             matrix of explanatory variables
  # coeffs        vector of coefficients, containing (in that order) coefficients for X and thresholds mu between classes
  # smooth_coeffs position of coefficients to be smoothed in the vector of coefficients
  # lambda        smoothing parameter
  # details       Boolean to indicate if the different components of the likelihood should be returned as well
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   LL          (negative) LL at given coefficients
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  k <- dim(X)[2]
  n <- length(coeffs)
  classes <- n-k+1
  mu <- c(-Inf,coeffs[(k+1):n],Inf)
  Y_star <- X %*% coeffs[1:k]
  p <- Y_star*0
  
  for (n in 1:classes) {
    i <- Y==n
    p[i] <- pnorm(mu[n+1]-Y_star[i])-pnorm(mu[n]-Y_star[i])
    
  }
  logp_coeffs <- sum(log(dnorm(sqrt(lambda) * diff(coeffs[smooth_coeffs]))))
  LL <- sum(log(p)) + logp_coeffs
  
  if (details){
    out <- NULL
    out$LL <- -LL
    out$LL_data <- sum(log(p))
    out$LL_smooth <- logp_coeffs
    out$p <- p
    return(out)
  }
  else {return(-LL)}
}