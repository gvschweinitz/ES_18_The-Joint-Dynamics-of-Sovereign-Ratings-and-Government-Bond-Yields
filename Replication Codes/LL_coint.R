LL_coint <- function(Ylin,Yoprob,Xlin,Xoprob,XEC,coeffs,NEC,Noprob,lambda,details=FALSE){
  # DESCRIPTION
  # This function calculates the joint log-likelihood of the linear and the oprob model given lambda
  # Estimation is performed without a call of polr, avoiding overhead costs due to multiple optimizer calls
  # The log-likelihood has the components:
  #     - fitting part, summed over observations (linear and oprob)
  #     - smoothness at (differenced) smoothed coefficients (linear and oprob, weighted)
  #----------------------------------------------------
  # USAGE
  # LL_coint(Ylin,Yoprob,Xlin,Xoprob,XEC,coeffs,NEC,Noprob,lambda,details=FALSE)
  #----------------------------------------------------
  #INPUT
  # Ylin          vector of the dependent linear variable
  # Yoprob        vector of the dependent oprob variable  (discrete number of values == classes)
  # Xlin          matrix of explanatory variables (linear)
  # Xoprob        matrix of explanatory variables (oprob)
  # XEC           matrix of longrun explanatory varialbes (equal relationship), yields being the first
  # coeffs        vector of coefficients to be optimized, containing long-run, short-run oprob, threshold oprob coefficients
  # NEC           number of variables in the cointegration matrix (less one)
  # Noprob        number of short-run oprob coefficients (without thresholds)
  # smooth_coeffs position of coefficients to be smoothed in the vector of coefficients
  # lambda        smoothing parameter
  # details       Boolean to indicate if the different components of the likelihood should be returned as well
  #----------------------------------------------------
  #OUTPUT
  # list with
  #   LL                  (negative) of overall log-likelihood
  #   LL_lin              log-likelihood of linear model
  #   LL_oprob            log-likelihood of ordered probit model
  #   LL_smooth_lin       log-likelihood-adjustment for smoothing in linear model
  #   LL_smooth_oprob     log-likelihood-adjustment for smoothing in ordered-probit model
  #   coefficients_lin    coefficients in linear equation
  #   coefficients_oprob  coefficients in ordered-probit equation
  #   coefficients_EC     coefficients of error-correction term
  #   lambda              smoothing coefficient
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  classes <- 3
  coeffs_EC <- coeffs[1:NEC]
  EC <- XEC %*% c(1,-coeffs_EC)
  Xest <- cbind(c(EC,sqrt(lambda)*diff(coeffs_EC)),Xlin)
  cest <- solve(t(Xest) %*% Xest) %*% (t(Xest) %*% Ylin)
  resids_lin <- Ylin - (Xest %*% cest)
  T <- dim(XEC)[1]
  sd_lin <- sd(resids_lin[1:T])
  LL_lin <- -T/2 * log(sd_lin^2)
  
  coeffs_oprob <- coeffs[(NEC+1):(NEC+Noprob)]
  Xest <- cbind(EC,Xoprob)
  mu <- c(-Inf,coeffs_oprob[(Noprob-1):Noprob],Inf)
  Y_star <- Xest %*% coeffs_oprob[1:(Noprob-2)]
  p <- Y_star*0
  for (state in 1:classes){
    n <- (Yoprob==(state-2))
    p[n] <- pnorm(mu[state+1]-Y_star[n])-pnorm(mu[state]-Y_star[n])
  }
  LL_oprob <- sum(log(p))
  
  # penalize differences for coefficients of facratings
  # penalization scaled down by regression alpha
  LL_smooth_lin <- sum(log(dnorm(cest[1] * (sqrt(lambda)/sd_lin) * diff(coeffs_EC))))
  LL_smooth_oprob <- sum(log(dnorm(coeffs_oprob[1] * sqrt(lambda) * diff(coeffs_EC))))
  LL <- LL_lin + LL_oprob + LL_smooth_lin + LL_smooth_oprob

  if (details){
    rownames(cest)[1] <- "EC"
    names(coeffs_oprob) <- c("EC",colnames(Xoprob))
    names(coeffs_EC) <- colnames(XEC)[-1]
    out <- NULL
    out$LL <- -LL
    out$LL_lin <- LL_lin
    out$LL_oprob <- LL_oprob
    out$LL_smooth_lin <- LL_smooth_lin
    out$LL_smooth_oprob <- LL_smooth_oprob
    out$coefficients_lin <- cest
    out$coefficients_oprob <- coeffs_oprob
    out$coefficients_EC <- c(1,-coeffs_EC)
    out$lambda <- lambda
    return(out)
  }
  else {
    return(-as.numeric(LL))
  }
}