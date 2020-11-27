compare_models <- function(y1,r1,y2=NULL,r2=NULL){
  # DESCRIPTION:
  # Performs model comparisons of two models. Can be one of:
  #   - separate models (y1, r1) versus cointegration model (y2)
  #   - separate models (y1, r1) versus separate models (y2, r2)
  #   - cointegration model (y1) versus cointegration model(r1)
  #-------------------------------------------------------------------------------
  #USAGE
  # compare_models(y1,r1,y2=NULL,r2=NULL)
  #-------------------------------------------------------------------------------
  #INPUT
  # y1    First (sub-)model
  # r1    Second (sub-)model
  # y2    Third (sub-)model, may be NULL
  # r2    Fourth (sub-)model, may be NULL
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   out$LLR           LLR of the test on the null_model
  #   out$k             DF of the test
  #   out$null_model    description of which of the models is put as the null models
  #   out$p             rejection probability of the null model
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  if (is.null(y2)){
    LLR <- 2*(y1$post$LL_out$LL_data-r1$post$LL_out$LL_data)
    k <- y1$k - r1$k
  }
  else{
    if (is.null(r2)){
      LLR <- 2*(y1$post$LL_out$LL_data+r1$post$LL_out$LL_data-y2$post$LL_out$LL_data)
      k <- y1$k + r1$k - y2$k
    }
    else{
      LLR <- 2*(y1$post$LL_out$LL_data+r1$post$LL_out$LL_data-y2$post$LL_out$LL_data-r2$post$LL_out$LL_data)
      k <- y1$k + r1$k - y2$k - r2$k
    }
  }
  
  null_model <- "second"
  if (k<0){
    LLR <- -LLR
    k <- -k
    null_model <- "first"
  }
  
  out <- NULL
  out$LLR <- LLR
  out$k <- k
  out$null_model <- null_model
  out$p <- 1-pchisq(LLR,k)
  
  return(out)
}