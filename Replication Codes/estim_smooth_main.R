estim_smooth_main <- function(pre,formula,lambda = 1,oprob = FALSE, calcgr = FALSE, calchess = FALSE, usemulti = FALSE){
  # DESCRIPTION:
  # Performs main estimation for a linear model or an ordered probit
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_smooth_main(pre,formula,lambda=1,oprob = FALSE, calcgr=FALSE, calchess=FALSE,usemulti=FALSE)
  #-------------------------------------------------------------------------------
  #INPUT
  # pre             result from pre-estimation (estimate_smooth_pre.R)
  # formula         formula to be estimated
  # oprob           boolean to indicate ordered probit
  # lambda          smoothing parameter for facratings (only needed to check if linear or full estimation is to be done)
  # calcgr          boolean to indicate if the final gradient should be calculated
  # calchess        boolean to indicate if the hessian should be calculated
  # usemulti        boolean if multiple optimization routines should be used
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   main$T        total number of observations
  #   main$N        total number of coefficients
  #   main$lambda   smoothing parameter
  #   main$kappa    effective number of coefficients in smoothing part
  #   main$k        effective number of coefficients
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
  
  lambdathresh <- 1000^3
  if (!is.formula(formula)) {
    formula <- as.formula(formula)
  }
  
  
  N <- length(pre$coeffs_start)
  T <- dim(pre$X)[1]
  Y <- pre$Y
  X <- pre$X
  
  base <- substr(pre$RHS,1,10)
  smooth_coeffs <- (base=="facratings")
  N_smooth <- sum(smooth_coeffs)
  X_smooth <- pre$X[,pre$RHS[smooth_coeffs]]
  D <- mat.or.vec(N_smooth-1,N_smooth)
  for (n in 1:(N_smooth-1)) {
    D[n,n] <- 1
    D[n,n+1] <- -1
  }
  D <- D*sqrt(lambda)
  
  # In order to avoid the creation of a TxT-matrix, get the trace directly
  Plambda_adj <- solve(t(X_smooth) %*% X_smooth + t(D) %*% D)
  kappa <- 0
  for (n in 1:T){
    kappa <- kappa + t(X_smooth[n,]) %*% Plambda_adj %*% X_smooth[n,]
  }
  k <- N-N_smooth+kappa
  
  if (oprob) {
    print("Calculating smoothed results")
    smooth_coeffs <- c(smooth_coeffs,FALSE,FALSE)
    coeffs <- pre$coeffs_start
    func <- function(x) {LL_oprob(Y = Y,X = X,coeffs = x,smooth_coeffs = smooth_coeffs,lambda = lambda)}
    gr_func <- function(x) {gr_oprob(Y = Y,X = X,coeffs = x,smooth_coeffs = smooth_coeffs,lambda = lambda)}
    hess_func <- function(x) {hess_oprob(Y = Y,X = X,coeffs = x,smooth_coeffs = smooth_coeffs,lambda = lambda)}
    if (lambda<lambdathresh){
      LL <- func(coeffs)
      if (!is.finite(LL)){
        print("setting coeffs for facratings to mean")
        coeffs[smooth_coeffs] <- mean(coeffs[smooth_coeffs],na.rm=TRUE)
        LL <- func(coeffs)
      }      
      if (!is.finite(LL)){
        print("setting coeffs to 0")
        coeffs <- coeffs*0
        coeffs[length(coeffs)-1] <- qnorm(mean(Y==1))
        coeffs[length(coeffs)] <- 1-qnorm(mean(Y==3))
        LL <- func(coeffs)
      }
      
      N_par <- length(coeffs)
      if (usemulti) {
        res_all = tryCatch({optimx(par = as.vector(coeffs),fn = func, gr=gr_func, 
                          method = c('Nelder-Mead', 'BFGS', 'nlm', 'nlminb', 'ucminf'))},
                           error = function(err){
                             print("error encountered - trying again without gradient")
                             print(err)
                             res <- optimx(par = as.vector(coeffs),fn = func, 
                                           method = c('Nelder-Mead', 'BFGS', 'nlm', 'nlminb', 'ucminf'))
                             return(res)}
        )
      }
      else{
        res_all = tryCatch({optimx(par = as.vector(coeffs),fn = func, gr=gr_func, 
                                   method = 'nlminb')},
                           error = function(err){
                                   print("error encountered - trying again without gradient")
                                   print(err)
                                   res <- optimx(par = as.vector(coeffs),fn = func, method = 'nlminb')
                                   return(res)}
        )
      }
      
      method <- rownames(res_all)[res_all[,"value"]==min(res_all[,"value"])]
      res_smooth <- NULL
      res_smooth$coefficients <- c(as.matrix(res_all[method,1:N_par]))
      res_smooth$fitted.values <-  pre$X %*% as.vector(res_smooth$coefficients[1:(N_par-2)])
      res_smooth$fitted.probs <- cbind(pnorm(res_smooth$coefficients[N_par-1]-res_smooth$fitted.values),
                                       1-pnorm(res_smooth$coefficients[N_par-1]-res_smooth$fitted.values)-pnorm(res_smooth$fitted.values-res_smooth$coefficients[N_par]),
                                       pnorm(res_smooth$fitted.values-res_smooth$coefficients[N_par]))
      names(res_smooth$coefficients) <- names(coeffs)
      coeffs <- c(res_smooth$coefficients)
    }
    else{
      res_smooth <- NULL
      res_smooth$coefficients <- pre$coeffs_start
      
      res_all <- NULL
      method <- NULL
    }
    
    if (calcgr){
      res_smooth$gr <- gr_func(coeffs)
      names(res_smooth$gr) <- names(coeffs)
    }
    else{res_smooth$gr<-NULL}
    if (calchess){
      print("Calculating Hessian")
      res_smooth$hess <- hess_func(coeffs)
      colnames(res_smooth$hess) <- names(coeffs)
      rownames(res_smooth$hess) <- names(coeffs)
    }
    else{res_smooth$hess <- NULL}
  }
  else{
    data_revised <- as.data.frame(cbind(pre$Y,pre$X))
    colnames(data_revised)[1] <- as.character(formula)[2]
    add_mat <- mat.or.vec(N_smooth-1,N+1)
    colnames(add_mat) <- colnames(data_revised)
    names <- pre$RHS[smooth_coeffs]
    for (n in 1:N_smooth) {
      add_mat[,names[n]] <- D[,n]
    }
    data_revised <- rbind(data_revised,add_mat)
    res_smooth <- lm(formula, data_revised)
    res_all <- NULL
    method <- NULL
  }
  
  main <- NULL
  main$T <- T
  main$N <- N
  main$lambda <- lambda
  main$kappa <- kappa
  main$k <- k
  main$smooth_coeffs <- smooth_coeffs
  main$res_smooth <- res_smooth
  main$res_all <- res_all
  main$method <- method
  return(main)
}