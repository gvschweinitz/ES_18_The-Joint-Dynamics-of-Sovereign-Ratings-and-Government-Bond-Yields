gr_oprob <- function(Y,X,coeffs,smooth_coeffs,lambda) {
  # DESCRIPTION
  # This function calculates the gradient log-likelihood for a ordered probit with (partially) smoothed coefficients.
  # The gradient has the components:
  #     - fitting part, summed over observations
  #     - fitting part, for thresholds
  #     - smoothness at (differenced) smoothed coefficients
  #----------------------------------------------------
  # USAGE
  # gr_oprob(Y,X,coeffs,smooth_coeffs,lambda)
  #----------------------------------------------------
  # INPUT
  # Y             vector of the dependent variable (discrete number of values == classes)
  # X             matrix of explanatory variables
  # coeffs        vector of coefficients, containing (in that order) coefficients for X and thresholds mu between classes
  # smooth_coeffs position of coefficients to be smoothed in the vector of coefficients
  # lambda        smoothing parameter
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   gr          Gradient of log-likelihood at given coefficients
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  t <- dim(X)[1]
  k <- dim(X)[2]
  cnt_coeffs <- length(coeffs)
  classes <- cnt_coeffs-k+1
  mu <- c(-Inf,coeffs[(k+1):cnt_coeffs],Inf)
  gr <- coeffs*0
  
  Y_star <- X %*% coeffs[1:k]
  Y_star_mu <- mat.or.vec(t,2)
  probs <- mat.or.vec(t,1)
  dens <- mat.or.vec(t,2)
  for (n in 1:classes){
    i <- (Y==n)
    if (n<classes){
      Y_star_mu[i,2] <- mu[n+1]-Y_star[i]
    }
    if (n>1) {
      Y_star_mu[i,1] <- mu[n]-Y_star[i]
    } 
  }
  dens <- dnorm(Y_star_mu)
  dens[Y_star_mu==0] <- 0
  p <- pnorm(Y_star_mu)
  p[Y_star_mu[,1]==0,1] <- 0
  p[Y_star_mu[,2]==0,2] <- 1
  probs <- p[,2]-p[,1]
  
  
  fac_beta <- mat.or.vec(t,1)
  fac_mu <- mat.or.vec(classes+1,1)
  for (n in 1:classes){
    i <- (Y==n)
    fac_beta[i] <- -(dens[i,2]-dens[i,1])/probs[i]
    fac_mu[n] <- fac_mu[n]-sum(dens[i,1]/probs[i])
    fac_mu[n+1] <- fac_mu[n+1]+sum(dens[i,2]/probs[i])
  }
  
  
  coeffs_s <- coeffs[smooth_coeffs]
  coeffs_s <- c(coeffs_s[1],coeffs_s,coeffs_s[length(coeffs_s)])
  
  gr[smooth_coeffs] <- lambda*diff(coeffs_s,differences=2)
  gr[1:k] <- gr[1:k] + fac_beta %*% X
  gr[(k+1):cnt_coeffs] <- fac_mu[2:classes]
  
  return(-gr)
}