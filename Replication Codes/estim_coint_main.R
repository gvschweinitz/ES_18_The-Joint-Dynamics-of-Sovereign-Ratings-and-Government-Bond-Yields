estim_coint_main <- function(pre,lambda = 1,usemulti = FALSE){
  # DESCRIPTION:
  # Performs main estimation for a cointegration model (without returning gradient and hessian)
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_coint_main(pre,formula,lambda=1,usemulti=FALSE)
  #-------------------------------------------------------------------------------
  #INPUT
  # pre             result from pre-estimation (estimate_coint_pre.R)
  # lambda          smoothing parameter for facratings (only needed to check if linear or full estimation is to be done)
  # usemulti        boolean if multiple optimization routines should be used
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   main$T        total number of observations
  #   main$N        total number of coefficients
  #   main$lambda   smoothing parameter
  #   main$kappa    reduction due to smoothing
  #   main$k        effective degrees of freedom
  #   main$smooth_coeffs  boolean vector giving places of coefficients to be smoothed
  #   main$res_smooth     list with result of the smoothed estimation. Can either be the result of an:
  #                         lm-estimation (call summary)
  #                         oprob estimation with entries "coefficients", "gr" and "hess"
  #   main$res_all  matrix output of optimx. Otherwise NULL
  #   main$method   optimal method from optimx. Otherwise NULL
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  library(matrixcalc)
  library(optimx)
  library(numDeriv)
  
  N <- dim(pre$Xlin)[2] + dim(pre$Xoprob)[2] + 2 + dim(pre$XEC)[2] + 1
  T <- dim(pre$Xlin)[1]
  
  eccols <- pre$RHSEC
  base <- substr(eccols,1,10)
  smooth_coeffs <- (base=="facratings")
  N_smooth <- sum(smooth_coeffs)
  X_smooth <- as.matrix(pre$XEC[,eccols[smooth_coeffs]])
  D <- mat.or.vec(N_smooth-1,N_smooth)
  for (n in 1:(N_smooth-1)) {
    D[n,n] <- 1
    D[n,n+1] <- -1
  }
  
  Plambda_adj <- solve(t(X_smooth) %*% X_smooth + lambda * t(D) %*% D)
  kappa <- 0
  for (n in 1:T){
    kappa <- kappa + t(X_smooth[n,]) %*% Plambda_adj %*% X_smooth[n,]
  }
  k <- N-N_smooth+kappa
  
  func <- function(x) {LL_coint(Ylin=Ylin,Yoprob=pre$Yoprob,Xlin=Xlin,Xoprob=pre$Xoprob,XEC=pre$XEC,coeffs=x,
                                 NEC=NEC,Noprob=Noprob,lambda=lambda,details=FALSE)}
  func_long <- function(x) {LL_coint(Ylin=Ylin,Yoprob=pre$Yoprob,Xlin=Xlin,Xoprob=pre$Xoprob,XEC=pre$XEC,coeffs=x,
                                     NEC=NEC,Noprob=Noprob,lambda=lambda,details=TRUE)}
  
  coeffs <- c(pre$coeffs_EC,pre$coeffs_oprob)
  NEC <- length(pre$coeffs_EC)
  Xlin <- as.matrix(rbind(pre$Xlin,mat.or.vec(NEC-1,dim(pre$Xlin)[2])))
  Ylin <- c(pre$Ylin,mat.or.vec(NEC-1,1))
  Noprob <- length(pre$coeffs_oprob)
  print("Calculating smoothed results")
  print(func(coeffs))
  if (usemulti) {res_all <- optimx(par = coeffs,fn = func,method = c('Nelder-Mead', 'BFGS', 'nlm', 'nlminb', 'ucminf')) }
  else {res_all <- optimx(par = coeffs,fn = func,method = c('nlminb'))}
  
  method <- rownames(res_all)[res_all[,"value"]==min(res_all[,"value"])]
  detailed_res <- func_long(as.matrix(res_all[method,1:(NEC+Noprob)]))
  print(paste("optimized at LL=",min(res_all[,"value"])))
  
  main <- NULL
  main$T <- T
  main$N <- N
  main$lambda <- lambda
  main$kappa <- kappa
  main$k <- k
  main$smooth_coeffs <- smooth_coeffs
  main$detailed_res <- detailed_res
  main$res_all <- res_all
  main$method <- method
  return(main)
}