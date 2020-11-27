bootstrap_par <- function(data,res_yields, res_ratings, eq_yields,eq_ratings,fspecs=fspecs, sample_info = NULL, NULLSTART = FALSE, iterations=1000) {
  #DESCRIPTION:
  # Estimates a set of "iteration" coefficients from bootstrapped data, using a wild bootstrap.
  #-------------------------------------------------------------------------------
  #USAGE
  # bootstrap_par(data,res_yields, res_ratings, eq_yields,eq_ratings,fspecs=fspecs, sample_info = NULL, NULLSTART = FALSE, iterations=1000)
  #-------------------------------------------------------------------------------
  #INPUT
  # data            data frame with data with the naming convention from data_prep
  # res_yields      results from yield estimation
  # res_ratings     results from rating estimation
  # eq_yields       formula for yield estimation
  # eq_ratings      formula for rating estimation
  # fspecs          Structure controlling estimation setup
  # sample_info     Result from analyze_sample.R
  # NULLSTART       Set starting coefficients for rating equation to zero, and constants to qnorm(1/3) and qnorm(2/3)
  # iterations      Number of iterations
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   tab_y         iterations x Ncoeffs_y matrix of bootstrap coefficients from yield equation
  #   tab_3         iterations x Ncoeffs_r matrix of bootstrap coefficients from rating equation
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz 
  library(dummies)
  library(foreach)
  library(doParallel)
  print("opening clusters")
  cl<-makeCluster(8, outfile="")
  registerDoParallel(cl)
  functionlist <- c("bootstrap_sample","simulation","extract_scenBS",
                    "estim_smooth","estim_smooth_pre","estim_smooth_main","estim_smooth_post","LL_oprob","gr_oprob","hess_oprob")
  
  # Set lag orders
  
  ylags_y = fspecs$ylags_y
  ylags_r = fspecs$ylags_r
  rlags_y = fspecs$rlags_y
  rlags_r = fspecs$rlags_r
  asym = fspecs$asym
  cdums <- fspecs$cdums
  rat_var = fspecs$idratings
  yield_var = fspecs$idyields_final
  crossid = fspecs$idc
  dateid = fspecs$idobs
  
  run_bs <- function(i){
    print(paste("run",i))
    check <- FALSE
    throw <- 1
    while (!check){
      new_sample <- bootstrap_sample(sample_info = sample_info, empirical_changes = empirical_changes,data=data,res_y=res_y,res_r=res_r,fspecs=fspecs)
      if (min(new_sample$data_r[,"lgratings2"]>min_lgratings2)){
        print(min(new_sample$data_r[,"lgratings2"]))
        throw <- throw+1
      }
      else {check <- TRUE}
    }
    res_y_bs <- estim_smooth(eq_yields,new_sample$data_y,lambda=lambda_y,oprob=FALSE, calcgr=FALSE, calchess=FALSE, usemulti=FALSE, polr_start = NULL, calcpost=FALSE)
    res_r_bs <- estim_smooth(eq_ratings,new_sample$data_r,lambda=lambda_r,oprob=TRUE, calcgr=FALSE, calchess=FALSE, usemulti=FALSE, polr_start=starter, calcpost=FALSE)
    out <- c(res_y_bs$res_smooth$coefficients,res_r_bs$res_smooth$coefficients,throw)
    return(out)
  }
  
  print("setting up bootstrap")
  l_thres = 1000^3
  blocks <- 10
  bw <- ceiling(iterations/blocks)
  t1 <- 1
  t2 <- bw
  
  # Setup wild bootstrap
  dist <- c(-sqrt(2),-1,-sqrt(0.5),sqrt(0.5),1,sqrt(1.5))
  N_dist = 6

  # Redo estimation
  if (!is.list(res_yields) | !is.list(res_ratings)){
    lambda_y <- res_yields
    lambda_r <- res_ratings
    res_y <- estim_smooth(eq_yields,data,lambda=lambda_y,oprob=FALSE, calcgr=FALSE, calchess=TRUE, usemulti=FALSE, polr_start = NULL, calcpost=FALSE)
    res_r <- estim_smooth(eq_ratings,data,lambda=lambda_r,oprob=TRUE, calcgr=FALSE, calchess=FALSE, usemulti=FALSE, polr_start=NULL, calcpost=FALSE)
  }
  else{
    res_y <- res_yields
    res_r <- res_ratings
    lambda_y <- res_yields$lambda
    lambda_r <- res_ratings$lambda
  }
  
  p_y <- length(res_y$res_smooth$coefficients)
  p_r <- length(res_r$res_smooth$coefficients)
  
  long <- mat.or.vec(iterations,p_y+p_r+1)
  colnames(long) <- c(names(res_y$res_smooth$coefficients),names(res_r$res_smooth$coefficients),"throw")
  
  # Extract residuals & legitimate selections starting scenario
  if (is.null(sample_info)) {sample_info <- analyze_sample(data,res_y, crossid = "country",dateid = "dateid",maxlag = max(rlags_r,rlags_y,ylags_r,ylags_y))}
  empirical_changes <- extract_changes(data)
  
  starter = res_r$res_smooth$coefficients
  print(starter)
  if (NULLSTART) {
    starter = starter * 0
    starter[length(starter)] = qnorm(2/3)
    starter[length(starter)-1] = qnorm(1/3)
  }
  
  min_lgratings2 <- 2
  
  print("running bootstrap")

  for (b in 1:blocks){
    print(paste("starting block",b,"(iterations",t1,"to",t2,")"))
#     long_b <- foreach (i=t1:t2,.combine=rbind,.export=functionlist) %do% (run_bs(i))
    long_b <- foreach (i=t1:t2,.combine=rbind,.export=functionlist) %dopar% (run_bs(i))
    long[t1:t2,] <- long_b
    t1<-t2+1
    t2<-min(iterations,t2+bw)
  }  
  
  stopCluster(cl)
  
  results <- NULL
  results$tab_y <- long[,1:p_y]
  results$tab_r <- long[,(p_y+1):(p_y+p_r)]
  throw <- sum(long[,"throw"])
  print(paste("total number of draws: ",throw,"(share rejections: ",round((throw-iterations)/iterations,2),")"))
  return(results)
}