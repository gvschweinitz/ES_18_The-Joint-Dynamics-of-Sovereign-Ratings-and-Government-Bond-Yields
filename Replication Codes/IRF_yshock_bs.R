IRF_yshock_bs <- function(data,tab_y,tab_r,std_y = 0,fspecs=fspecs,empirical_changes=NULL,start=NULL,scen=NULL,periods=120,diff_to_benchmark=TRUE){
  #DESCRIPTION:
  # Calculates IRFs based on bootstrap coefficient draws for a given scenario and given initial values for yields and ratings.
  # For every bootstrap draw, the median IRF is calculated taking upgrade and downgrade uncertainty into account (run simulation.R 1000 times)
  #-------------------------------------------------------------------------------
  #USAGE
  # IRF_yshock_bs(data,tab_y,tab_r,std_y = 0,fspecs=fspecs,empirical_changes=NULL,start=NULL,scen=NULL,periods=120,diff_to_benchmark=TRUE)
  #-------------------------------------------------------------------------------
  #INPUT
  # data              data frame with data with the naming convention from data_prep
  # tab_y             bootstrap coefficients from yield equation
  # tab_r             bootstrap coefficients from rating equation
  # std_y             Standard deviation of yield changes. Needed to account for the influence of past changes due to asymmetric effects
  # fspecs            Structure controlling estimation setup
  # empirical_changes Result from extract_changes.R
  # start             Matrix with columns "ratings" and "yields". Can be used to construct scenario scen
  # scen              IRF scenario. If NULL, scenario is constructed from start using "build_scen.R" as the changes since start
  # periods           IRF horizon
  # diff_to_benchmark Boolean. FALSE: IRF relative to initial rating. TRUE: IRF relative to usual adjustment towards long-run equilibrium
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   IRF_y           ndraws x (periods+length(downgrade)+1) matrix of IRFs of yields
  #   IRF_r           ndraws x (periods+length(downgrade)+1) matrix of IRFs of ratings
  #   Default         ndraws x (periods+length(downgrade)+1) binary matrix of IRFs of defaults
  #   scen            IRF scenario. If NULL, scenario is constructed from start using "build_scen.R" as the changes since start
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  library(foreach)
  library(doParallel)
  print("opening clusters")
  cl<-makeCluster(8, outfile="")
  registerDoParallel(cl)
  functionlist <- c("simulation")
  
  minrat <- min(fspecs$factors)-1
  if (is.null(start)){
    yields <- c(20,20)
    ratings <- c(9,7)
  }
  else {
    yields <- start$yields
    ratings <- start$ratings
  }
  
  if (!diff_to_benchmark){
    table_bm <- mat.or.vec(periods,2)
    colnames(table_bm) <- c("Yields","Rating")
  }

  #number of periods in the IRF
  res_periods = periods+length(yields)
  # Build scenario based on initial yield and rating scenario
  if (is.null(scen)){scen <- build_scen(yields,ratings,tab_y[1,],tab_r[1,],fspecs = fspecs)}
  # Build scenario based on first yield and rating in the scenario, equal length as scen
  if (diff_to_benchmark){scen_benchmark <- build_scen(yields=rep(yields[1],length(yields)),ratings=rep(ratings[1],length(ratings)),tab_y[1,],tab_r[1,],fspecs = fspecs)}
  
  # number of repetitions needed to get stable median
  reps_sim <- 1000
  reps <- dim(tab_y)[1]
  IRF_y <- mat.or.vec(reps,res_periods)
  IRF_r <- mat.or.vec(reps,res_periods)
  Default <- mat.or.vec(reps,res_periods)

  calc_IRF <- function(i){
    # Function for calculating the IRF for one coefficient vector
    # Calculates it as the median IRF with respect to rating shocks and possible yield shocks
    # In case of diff_to_benchmark = TRUE, differences to the benchmark are applied before taking the median
    coefficients_y <- tab_y[i,]
    coefficients_r <- tab_r[i,]
    y_sim <- mat.or.vec(reps_sim,res_periods)
    r_sim <- mat.or.vec(reps_sim,res_periods)
    r_sim_lvl <- mat.or.vec(reps_sim,res_periods)
    if (std_y>0){y_shock = rnorm(periods,0,std_y)}
    else{y_shock=0}
    
    for (r in 1:reps_sim){
      table_scen <- simulation(coefficients_y, coefficients_r,fspecs=fspecs, start = scen, y_shock = 0, empirical_changes = NULL, periods = periods, output = NULL, details = FALSE, countries = NA)
      if (diff_to_benchmark){
        table_bm <- simulation(coefficients_y, coefficients_r,fspecs=fspecs, start = scen_benchmark, y_shock = 0, empirical_changes = NULL, periods = periods, output = NULL, details = FALSE, countries = NA)
      }
      y_sim[r,] <- table_scen[,"Yield"] - table_bm[,"Yield"]
      r_sim[r,] <- table_scen[,"Rating"] - table_bm[,"Rating"]
      r_sim_lvl[r,] <- table_scen[,"Rating"]
    }
    ymed <- apply(y_sim,2,quantile,0.5)
    rmed <- apply(r_sim,2,quantile,0.5)
    isdefault <- apply(r_sim_lvl,2,quantile,0.5)==minrat
    return(rbind(ymed,rmed,isdefault))
  }
    
  print(paste("running bootstrap, calculating median using",reps_sim,"repetitions"))
  blocks <- 5
  bw <- ceiling(reps/blocks)
  t1 <- 1
  t2 <- bw
  for (b in 1:blocks){
    print(paste("starting block",b,"(reps",t1,"to",t2,")"))
    long_b <- foreach (i=t1:t2,.combine=rbind,.export=functionlist) %dopar% (calc_IRF(i))
    if (b==1){long<-long_b}
    else{long<-rbind(long,long_b)}
    t1<-t2+1
    t2<-min(reps,t2+bw)
  }  
  IRF_y <- long[seq(1,dim(long)[1]-2,3),]
  IRF_r <- long[seq(2,dim(long)[1]-1,3),]
  Default <- long[seq(3,dim(long)[1],3),]
  
  stopCluster(cl)
  results <- NULL
  results$IRF_y <- IRF_y
  results$IRF_r <- IRF_r
  results$Default <- Default
  results$scen <- scen
  return(results)
}