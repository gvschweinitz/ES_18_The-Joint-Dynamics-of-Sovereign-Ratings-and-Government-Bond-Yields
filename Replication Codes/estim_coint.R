estim_coint <- function(lambda=1,data,lags_yields,lags_ratings,contemp=TRUE,factors=c(7:24),constlong=FALSE,cdums=FALSE,rasym=TRUE,yasym=TRUE, calcpost=TRUE){
  # DESCRIPTION:
  # estimates a cointegrated model jointly
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_coint(lambda=1,data,lags_yields,lags_ratings,contemp=TRUE,factors=c(7:24),constlong=FALSE,cdums=FALSE,rasym=TRUE,yasym=TRUE, calcpost=TRUE)
  #-------------------------------------------------------------------------------
  #ARGUMENTS
  # lambda          smoothing parameter for facratings
  # data            data frame with the naming convention from data_prep
  # lags_yields     number of lags in the yields equation
  # lags_ratings    number of lags in the ratings equation
  # contemp         estimation with or without contemporaneous differences (Pesaran style)
  # factors         vector containing the rating dummies to be included
  # constlong       boolean to indicate if the constant should be included in the long-run relation
  # cdums           boolean if country dummies should be used
  # rasym           boolean if asymmetric effects are allowed in the rating equation
  # yasym           boolean if asymmetric effects are allowed in the yield equation
  # calcpost        boolean if post-estimation analytics should be performed
  #-------------------------------------------------------------------------------
  #OUTPUT
  # res is a list with entries
  #   pre           result of pre-estimation
  #   post          result of post-estimation
  #   several entries from main estimation
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  pre <- estim_coint_pre(data,fspecs=fspecs,constlong = constlong)
  res <- estim_coint_main(pre, lambda=lambda, usemulti=FALSE)
  if (calcpost){
    res$post <- estim_coint_post(pre,res)
  }
  res$pre <- pre
  
  return(res)
}