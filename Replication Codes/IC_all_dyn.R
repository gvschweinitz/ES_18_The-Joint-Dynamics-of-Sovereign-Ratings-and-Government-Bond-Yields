IC_all_dyn <- function(formula,data,max_diff = 1, oprob = FALSE, func = NULL, optimcol = 3) {
  # DESCRIPTION:
  # finds the (approximately) optimal smoothing parameter
  #-------------------------------------------------------------------------------
  #USAGE
  # IC_all_dyn(formula,data,max_diff = 1, oprob = FALSE, func = NULL, optimcol = 3)
  #-------------------------------------------------------------------------------
  #INPUT
  # formula         formula to be estimated
  # data            data frame with the naming convention from data_prep
  # max_diff        maximum difference between the optimal and the second-best lambda (stopping criterion)
  # oprob           boolean to indicate ordered probit
  # func            function for the calculation of initial lambda values
  # optimcol        column in the output table after which one should optimize
  #-------------------------------------------------------------------------------
  #OUTPUT
  # out             K x 5 matrix with columns lambda, AICc, BIC, kappa and LL_data
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  if (is.null(func)) {
    func <- function(x){exp(x)}
  }
  
  perform_est <- function(pre,formula,data,lambda,oprob){
    if(lambda>=lambda_threshold){
      pre <- estim_smooth_pre(formula=formula,data=data,oprob=oprob,lambda=lambda,polr_start=NULL)
    }
    main <- estim_smooth_main(pre=pre,formula=formula,lambda=lambda,oprob = oprob, calcgr=FALSE, calchess=FALSE,usemulti=FALSE)
    post <- estim_smooth_post(pre,main,oprob=oprob)
    res_vec <- c(lambda,post$AICc,post$BIC,main$kappa,post$LL_out$LL_data)
    return(res_vec)
  }
  
  lambda_threshold <- 1000^3
  max_l_start <- 8
  ncols <- 5
  table <- mat.or.vec(max_l_start,ncols)
  pre <- estim_smooth_pre(formula=formula,data=data,oprob=oprob,lambda=1,polr_start=NULL)
  for (l in 1:max_l_start) {
    print(l)
    lambda <- func(l)
    table[l,] <- perform_est(pre,formula,data,lambda,oprob)
  }
  colnames(table) <- c("lambda","AICc","BIC","kappa","LL_data")
  print(table)
  bestrow <- which(table[,optimcol]==min(table[,optimcol]))
  
  check <- (bestrow == 1)
  while (check){
    nrows <- dim(table)[1]
    if (table[1,"lambda"]==0){lambda <- 0.5*table[2,"lambda"]}
    else{lambda <- 0.5*table[1,"lambda"]}
    print(lambda)
    res_vec <- perform_est(pre,formula,data,lambda,oprob)
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
    lambda <- min(2*table[bestrow,"lambda"],lambda_threshold)
    print(lambda)
    res_vec <- perform_est(pre,formula,data,lambda,oprob)
    table <- rbind(table,res_vec)
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    check <- (bestrow == dim(table)[1])
    if (lambda>=lambda_threshold){
      if (check){return(table)}
    }
  }
  
  dt <- diff(table[,1])
  while (max(dt[(bestrow-1):min(bestrow,length(dt))])>max_diff){
    nrows <- dim(table)[1]
    lambda <- c(mean(table[(bestrow-1):bestrow,1]),mean(table[bestrow:(bestrow+1),1]))
    print(lambda)
    res_vec1 <- perform_est(pre,formula,data,lambda[1],oprob)
    res_vec2 <- perform_est(pre,formula,data,lambda[2],oprob)
    
    table <- rbind(table[1:(bestrow-1),],res_vec1,table[bestrow,],res_vec2,table[(bestrow+1):nrows,])
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    dt <- diff(table[,1])
  }
  
  return(table)
}