hess_oprob <- function(Y,X,coeffs,smooth_coeffs,lambda) {
  # DESCRIPTION
  # This function calculates the hessian of the log-likelihood for a ordered probit with (partially) smoothed coefficients.
  # The hessian has the components:
  #     - fitting part, summed over observations
  #     - fitting part, for thresholds
  #     - smoothness at (differenced) smoothed coefficients
  #----------------------------------------------------
  # USAGE
  # hess_oprob(Y,X,coeffs,smooth_coeffs,lambda)
  #----------------------------------------------------
  # ARGUMENTS
  # Y             vector of the dependent variable (discrete number of values == classes)
  # X             matrix of explanatory variables
  # coeffs        vector of coefficients, containing (in that order) coefficients for X and thresholds mu between classes
  # smooth_coeffs position of coefficients to be smoothed in the vector of coefficients
  # lambda        smoothing parameter
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   hess        Hessian of log-likelihood at given coefficients
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  t <- dim(X)[1]
  k <- dim(X)[2]
  cnt_coeffs <- length(coeffs)
  classes <- cnt_coeffs-k+1
  mu <- c(-Inf,coeffs[(k+1):cnt_coeffs],Inf)
  
  Y_star <- X %*% coeffs[1:k]
  Y_star_mu <- mat.or.vec(t,2)
  
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
  fac_betamu <- mat.or.vec(t,classes+1)
  hess_mu <- mat.or.vec(classes+1,classes+1)
  hess <- mat.or.vec(cnt_coeffs,cnt_coeffs)
  for (n in 1:classes){
    i <- (Y==n)
    fac_beta[i] <- (Y_star_mu[i,1]*dens[i,1]-Y_star_mu[i,2]*dens[i,2])/probs[i] - ((dens[i,2]-dens[i,1])/probs[i])^2
    fac_betamu[i,n] <- fac_betamu[i,n] - Y_star_mu[i,1]*dens[i,1]/probs[i] - (dens[i,1]*(dens[i,2]-dens[i,1]))/(probs[i]^2)
    fac_betamu[i,n+1] <- fac_betamu[i,n+1] + Y_star_mu[i,2]*dens[i,2]/probs[i] + (dens[i,2]*(dens[i,2]-dens[i,1]))/(probs[i]^2)
    
    hess_mu[n,n] <- hess_mu[n,n] + sum(Y_star_mu[i,1]*dens[i,1]/probs[i] - (dens[i,1]/probs[i])^2)
    hess_mu[n+1,n+1] <- hess_mu[n+1,n+1] - sum(Y_star_mu[i,2]*dens[i,2]/probs[i] + (dens[i,2]/probs[i])^2)
    hess_mu[n,n+1] <- hess_mu[n,n+1] + sum(dens[i,1]*dens[i,2]/(probs[i])^2)
    hess_mu[n+1,n] <- hess_mu[n,n+1]
  }
  
  for (i in 1:t){
    hess[1:k,1:k] <- hess[1:k,1:k]+fac_beta[i] * (X[i,]%*%t(X[i,]))
  }
  hess[(k+1):cnt_coeffs,1:k] <- t(fac_betamu[,2:classes]) %*% X
  hess[1:k,(k+1):cnt_coeffs] <- t(hess[(k+1):cnt_coeffs,1:k])
  hess[(k+1):cnt_coeffs,(k+1):cnt_coeffs] <- hess_mu[2:classes,2:classes]
  
  pos <- (1:cnt_coeffs)[smooth_coeffs]
  for (i in pos){
    hess[i,i] <- hess[i,i]-lambda*ifelse(i>pos[1]&i<pos[length(pos)],2,1)
    if (i>pos[1]){
      
      hess[i-1,i] <- hess[i-1,i]+lambda
      hess[i,i-1] <- hess[i-1,i]
    }
  }
  
  rownames(hess) <- rownames(coeffs)
  colnames(hess) <- rownames(coeffs)
  
  return(-hess)
}