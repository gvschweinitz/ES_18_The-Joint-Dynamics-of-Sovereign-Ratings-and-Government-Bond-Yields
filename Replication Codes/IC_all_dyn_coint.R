IC_all_dyn_coint <- function(data,fspecs=fspecs,constlong = FALSE,max_diff = 1, func = NULL, optimcol = 3) {
  # DESCRIPTION:
  # finds the (approximately) optimal smoothing parameter that holds similarly for the yield and the rating equation
  # Estimates a single long-run relationship in a cointegration framework
  #-------------------------------------------------------------------------------
  #USAGE
  # coint_BIC <- IC_all_dyn_coint(data,optlag_y,optlag_r,cdums=TRUE,rasym=TRUE,yasym=TRUE)
  #-------------------------------------------------------------------------------
  #INPUT
  # data            data frame with the naming convention from data_prep
  # fspecs          Structure controlling estimation setup
  # constlong       boolean to indicate if the constant should be included in the long-run relation
  # max_diff        maximum difference between the optimal and the second-best lambda
  # func            function for the calculation of initial lambda values
  # optimcol        column in the output table after which one should optimize
  #-------------------------------------------------------------------------------
  #OUTPUT
  # out             K x 7 matrix with columns lambda, AICc, BIC, kappa and LL_full, LL_lin and LL_oprob
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  if (is.null(func)) {
    func <- function(x){exp(x)}
  }
  
  perform_est <- function(pre,lambda){
    main <- estim_coint_main(pre=pre,lambda=lambda,usemulti=FALSE)
    post <- estim_coint_post(pre,main)
    res_vec <- c(lambda,post$AICc,post$BIC,main$kappa,main$detailed_res$LL,main$detailed_res$LL_lin,main$detailed_res$LL_oprob)
    return(res_vec)
  }
  
#   lambda_threshold <- 1000^3
  max_l_start <- 6
  ncols <- 7
  table <- mat.or.vec(max_l_start,ncols)
  pre <- estim_coint_pre(data,fspecs=fspecs,constlong = constlong)
  for (l in 1:max_l_start) {
    print(l)
    lambda <- func(l)
    table[l,] <- perform_est(pre,lambda)
  }
  colnames(table) <- c("lambda","AICc","BIC","kappa","LL_full","LL_lin","LL_oprob")
  print(table)
  bestrow <- which(table[,optimcol]==min(table[,optimcol]))
  
  check <- (bestrow == 1)
  while (check){
    nrows <- dim(table)[1]
    if (table[1,"lambda"]==0){lambda <- 0.5*table[2,"lambda"]}
    else{lambda <- 0.5*table[1,"lambda"]}
    print(lambda)
    res_vec <- perform_est(pre,lambda)
    if (table[1,"lambda"]==0){table <- rbind(table[1,],res_vec,table[2:nrows,])}
    else{table <- rbind(res_vec,table)}
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    check <- (bestrow == 1)
    
    if (table[2,"lambda"]-table[1,"lambda"]<max_diff){
      if(check){return(table)}
    }
  }

  check <- (bestrow == dim(table)[1])
  while (check){
    lambda <- 2*table[bestrow,"lambda"]
    print(lambda)
    res_vec <- perform_est(pre,lambda)
    table <- rbind(table,res_vec)
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    check <- (bestrow == dim(table)[1])
    if (lambda>=lambda_threshold){
      if (check){return(table)}
    }
  }
  
  print(table)
  dt <- diff(table[,1])
  while (max(dt[(bestrow-1):bestrow])>max_diff){
    nrows <- dim(table)[1]
    lambda <- c(mean(table[(bestrow-1):bestrow,1]),mean(table[bestrow:(bestrow+1),1]))
    print(lambda)
    res_vec1 <- perform_est(pre,lambda[1])
    res_vec2 <- perform_est(pre,lambda[2])
    
    table <- rbind(table[1:(bestrow-1),],res_vec1,table[bestrow,],res_vec2,table[(bestrow+1):nrows,])
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    dt <- diff(table[,1])
  }
  
  return(table)
}