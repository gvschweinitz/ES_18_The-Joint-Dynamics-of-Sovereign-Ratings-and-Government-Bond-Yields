IC_all_dyn_joint <- function(formula_y,formula_r,data,max_diff = 1, func = NULL, optimcol = 4) {
  # DESCRIPTION:
  # finds the (approximately) optimal smoothing parameter that holds similarly for the yield and the rating equation
  # However, it still estimates two separate long-run relationships!
  #-------------------------------------------------------------------------------
  #USAGE
  # IC_all_dyn_joint(formula,data,max_diff = 1, oprob = FALSE, func = NULL, optimcol = 3)
  #-------------------------------------------------------------------------------
  #INPUT
  # formula_y       formula of yield equation
  # formula_r       formula of rating equation
  # data            data frame with the naming convention from data_prep
  # max_diff        maximum difference between the optimal and the second-best lambda (stopping criterion)
  # func            function for the calculation of initial lambda values
  # optimcol        column in the output table after which one should optimize
  #-------------------------------------------------------------------------------
  #OUTPUT
  # out             K x 7 matrix with columns lambda, BICy, BICr, BICjoint, LL_data_y, LL_data_r, LL_joint
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  if (is.null(func)) {
    func <- function(x){exp(x)}
  }
  
  perform_est <- function(pre_y,pre_r,formula_y,formula_r,data,lambda){
    if(lambda>=lambda_threshold){
      pre_r <- estim_smooth_pre(formula=formula_r,data=data,oprob=TRUE,lambda=lambda,polr_start=NULL)
    }
    main_y <- estim_smooth_main(pre=pre_y,formula=formula_y,lambda=lambda,oprob = FALSE, calcgr=FALSE, calchess=FALSE,usemulti=FALSE)
    post_y <- estim_smooth_post(pre_y,main_y,oprob=FALSE)
    main_r <- estim_smooth_main(pre=pre_r,formula=formula_r,lambda=lambda,oprob = TRUE, calcgr=FALSE, calchess=FALSE,usemulti=FALSE)
    post_r <- estim_smooth_post(pre_r,main_r,oprob=TRUE)
    res_vec <- c(lambda,
                 post_y$BIC,
                 post_r$BIC,
                 post_y$BIC+post_r$BIC,
                 post_y$LL_out$LL_data,
                 post_r$LL_out$LL_data,
                 post_y$LL_out$LL_data+post_r$LL_out$LL_data)
    return(res_vec)
  }
  
  lambda_threshold <- 1000^3
  max_l_start <- 8
  ncols <- 7
  table <- mat.or.vec(max_l_start,ncols)
  colnames(table) <- c("lambda","BIC_y","BIC_r","BIC_joint","LL_data_y","LL_data_r","LL_joint")
  pre_y <- estim_smooth_pre(formula=formula_y,data=data,oprob=FALSE)
  pre_r <- estim_smooth_pre(formula=formula_r,data=data,oprob=TRUE)
  for (l in 1:max_l_start) {
    print(l)
    lambda <- func(l)
    res_vec <- perform_est(pre_y,pre_r,formula_y,formula_r,data,lambda)
    table[l,] <- res_vec
  }
  print(table)
  bestrow <- which(table[,optimcol]==min(table[,optimcol]))
  
  
  check <- (bestrow == 1)
  while (check){
    nrows <- dim(table)[1]
    if (table[1,"lambda"]==0){lambda <- 0.5*table[2,"lambda"]}
    else{lambda <- 0.5*table[1,"lambda"]}
    print(lambda)
    res_vec <- perform_est(pre_y,pre_r,formula_y,formula_r,data,lambda)
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
    res_vec <- perform_est(pre_y,pre_r,formula_y,formula_r,data,lambda)
    table <- rbind(table,res_vec)
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    check <- (bestrow == dim(table)[1])
    if (lambda>=lambda_threshold){
      if (check){return(table)}
    }
  }
  
  print(table)
  dt <- diff(table[,1])
  while (max(dt[(bestrow-1):min(bestrow,length(dt))])>max_diff){
    nrows <- dim(table)[1]
    lambda <- c(mean(table[(bestrow-1):bestrow,1]),mean(table[bestrow:(bestrow+1),1]))
    print(lambda)
    res_vec1 <- perform_est(pre_y,pre_r,formula_y,formula_r,data,lambda[1])
    res_vec2 <- perform_est(pre_y,pre_r,formula_y,formula_r,data,lambda[2])
    
    table <- rbind(table[1:(bestrow-1),],res_vec1,table[bestrow,],res_vec2,table[(bestrow+1):nrows,])
    bestrow <- which(table[,optimcol]==min(table[,optimcol]))
    dt <- diff(table[,1])
  }
  
  return(table)
}