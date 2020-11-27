estim_smooth <- function(formula,data,lambda = 1,oprob = FALSE, calcgr=TRUE, calchess = TRUE, usemulti = FALSE, polr_start = NULL, calcpost=TRUE){
  # DESCRIPTION:
  # estimates either a linear model or an ordered probit
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_smooth(formula,data,lambda = 1,oprob = FALSE, calcgr=TRUE, calchess = TRUE, usemulti = FALSE, polr_start = NULL, calcpost=TRUE)
  #-------------------------------------------------------------------------------
  #ARGUMENTS
  # formula         formula to be estimated
  # data            data frame with the naming convention from data_prep
  # lambda          smoothing parameter for facratings
  # oprob           boolean to indicate ordered probit
  # calcgr          boolean to indicate if the gradient should be calculated
  # calchess        boolean to indicate if the hessian should be calculated
  # usemulti        boolean if multiple optimization routines should be used
  # polr_start      starting values for polr estimation (not necessary)
  # calcpost        boolean if post-estimation analytics should be performed
  #-------------------------------------------------------------------------------
  #OUTPUT
  # res is a list with entries
  #   pre           result of pre-estimation
  #   post          result of post-estimation
  #   several entries from main estimation
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  pre <- estim_smooth_pre(formula, data, oprob=oprob, lambda=lambda, polr_start=polr_start)
  res <- estim_smooth_main(pre, formula, lambda=lambda, oprob=oprob, calcgr=calcgr, calchess=calchess, usemulti=usemulti)
  if (calcpost){
    res$post <- estim_smooth_post(pre,res,oprob=oprob)
  }
  res$pre <- pre
  
  return(res)
}