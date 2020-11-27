IRF_scenario_par <- function(month,year,country,data,res_yields,res_ratings,coeffs_bs,fspecs=fspecs,periods_start=6,empirical_changes=NULL,conf_lvls=NULL){
  #DESCRIPTION:
  # Calculates IRFs based on bootstrap coefficient draws for a given real-world scenario
  # Compares the simulated development of ratings and yields with up to 120 months of observed developments
  # Accounts for normal shocks to yields in simulation
  #-------------------------------------------------------------------------------
  #USAGE
  # IRF_scenario_par(month,year,country,data,res_yields,res_ratings,coeffs_bs,fspecs=fspecs,periods_start=6,empirical_changes=NULL,conf_lvls=NULL){
  #-------------------------------------------------------------------------------
  #INPUT
  # month             Month of the first observation of the real-world scenario
  # year              Year of the first observation of the real-world scenario
  # country           Country of the real-world scenario
  # data              data frame with data with the naming convention from data_prep
  # res_yields        results from yield estimation
  # res_ratings       results from rating estimation
  # coeffs_bs         Bootstrap coefficients for both equations. Output of "bootstrap_par.R"
  # fspecs            Structure controlling estimation setup
  # periods_start     Length of the real-world scenario before simulation
  # empirical_changes Result from extract_changes.R
  # conf_lvls         Confidence levels for IRFs. May also contain strings indicating standard deviations -> median +/- one standard deviation
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   scen            IRF scenario
  #   IRF_med_y       1 x periods vector of median IRFs of yields
  #   IRF_med_r       1 x periods vector of median IRFs of ratings
  #   IRF_conf_med_y  periods x length(conf_lvl) matrix of IRFs of yields at all conf_lvls
  #   IRF_conf_med_r  periods x length(conf_lvl) matrix of IRFs of ratings at all conf_lvls
  #   conf_lvls       Confidence levels for IRFs.  
  #   month           Month of the first observation of the real-world scenario
  #   year            Year of the first observation of the real-world scenario
  #   country         Country of the real-world scenario
  #   periods_start   Length of the real-world scenario before simulation
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz

  library(foreach)
  library(doParallel)
  print("opening clusters")
  cl<-makeCluster(8, outfile="")
  registerDoParallel(cl)
  functionlist <- c("simulation")
  
  if (is.null(conf_lvls)){conf_lvls <- c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99)}
  
  std_y <- sd(res_yields$res_smooth$residuals[1:length(res_yields$pre$Y)])
  
  tab_y <- coeffs_bs$tab_y
  tab_r <- coeffs_bs$tab_r
  reps_sim <- 1000
  reps <- dim(tab_y)[1]
  rat_var <- fspecs$idratings
  yield_var <- fspecs$idyields_final
  # ------------------------
  print(paste("Starting scenario for",country,"in",month,"-",year))
  scen <- extract_scenario(month,year,country,periods=periods_start,data,res_yields,res_ratings,rat_var=rat_var,yield_var=yield_var)
  par(mfrow = c(1,2))
  plot(scen$obs_yield)
  abline(v=periods_start)
  plot(scen$obs_ratings)
  abline(v=periods_start)
  print(cbind(scen$yield,scen$ratings))
  #number of periods in the IRF
  periods <- scen$period_sim
  periods_sim <- periods-periods_start
  
  calc_IRF <- function(i){
    # Function for calculating the IRF for one coefficient vector
    # Calculates it as the median IRF with respect to rating shocks and possible yield shocks
    # In case of diff_to_benchmark = TRUE, differences to the benchmark are applied before taking the median
    coefficients_y <- tab_y[i,]
    coefficients_r <- tab_r[i,]
    y_sim <- mat.or.vec(reps_sim,periods)
    r_sim <- mat.or.vec(reps_sim,periods)
    
    for (r in 1:reps_sim){
      if (std_y>0){y_shock = rnorm(periods_sim,0,std_y)}
      else{y_shock=0}
      table_scen <- simulation(coefficients_y,coefficients_r,fspecs=fspecs,start = scen,y_shock = y_shock,empirical_changes = empirical_changes,periods = periods_sim,output = NULL, details = FALSE)
      y_sim[r,] <- table_scen[,"Yield"]
      r_sim[r,] <- table_scen[,"Rating"]
    }
    res <- NULL
    y_med <- apply(y_sim,2,quantile,0.5)
    r_med <- apply(r_sim,2,quantile,0.5)
    res <- rbind(y_sim,r_sim,med_y=y_med,med_r=r_med)
    return(res)
  }
  
  long_b <- foreach (i=1:reps,.combine=rbind,.export=functionlist) %dopar% (calc_IRF(i))
# long_b <- foreach (i=1:reps,.combine=rbind,.export=functionlist) %do% (calc_IRF(i))
  pos_y <- which(rownames(long_b)=="med_y")
  pos_r <- which(rownames(long_b)=="med_r")
  temp_med_y <- long_b[pos_y,]
  temp_med_r <- long_b[pos_r,]
  long_b <- long_b[-c(pos_y,pos_r),]
  pos_y <- mat.or.vec(reps*reps_sim,1)
  for (i in 1:reps){pos_y[((i-1)*reps_sim+1):(i*reps_sim)] <- ((i-1)*reps_sim+1):(i*reps_sim)+(i-1)*reps_sim}
  temp_all_y  <-long_b[pos_y,]
  temp_all_r  <-long_b[pos_y+reps_sim,]
  rm(long_b)
  
  print("Calculating median and confidence bands")
  med_y <- apply(temp_med_y,2,quantile,0.5)
  utility <- rep(0,dim(temp_med_y)[1])
  for (i in 1:dim(temp_med_y)[1]){
    utility[i] <- sum((temp_med_y[i,]-med_y)^2)
  }
  IRF_med_y <- temp_med_y[which(utility==min(utility))[1],]
  
  med_r <- apply(temp_med_r,2,quantile,0.5)
  utility <- rep(0,dim(temp_med_r)[1])
  for (i in 1:dim(temp_med_r)[1]){
    utility[i] <- sum((temp_med_r[i,]-med_r)^2)
  }
  IRF_med_r <- temp_med_r[which(utility==min(utility))[1],]
  
  IRF_conf_med_y <- mat.or.vec(periods,length(conf_lvls))
  IRF_conf_med_r <- mat.or.vec(periods,length(conf_lvls))
  IRF_conf_shock_y <- mat.or.vec(periods,length(conf_lvls))
  IRF_conf_shock_r <- mat.or.vec(periods,length(conf_lvls))
  colnames(IRF_conf_med_y) <- conf_lvls
  colnames(IRF_conf_med_r) <- conf_lvls
  colnames(IRF_conf_shock_y) <- conf_lvls
  colnames(IRF_conf_shock_r) <- conf_lvls
  
  for (lvl in 1: length(conf_lvls)){
    IRF_conf_med_y[,lvl] <- apply(temp_med_y,2,quantile,as.numeric(conf_lvls[lvl]))
    IRF_conf_med_r[,lvl] <- apply(temp_med_r,2,quantile,as.numeric(conf_lvls[lvl]))
    IRF_conf_shock_y[,lvl] <- apply(temp_all_y,2,quantile,as.numeric(conf_lvls[lvl]))
    IRF_conf_shock_r[,lvl] <- apply(temp_all_r,2,quantile,as.numeric(conf_lvls[lvl]))
  }

  res <- NULL
  res$scen <- scen
  res$IRF_med_y <- IRF_med_y
  res$IRF_med_r <- IRF_med_r
  res$IRF_conf_med_y <- IRF_conf_med_y
  res$IRF_conf_med_r <- IRF_conf_med_r
  res$IRF_conf_shock_y <- IRF_conf_shock_y
  res$IRF_conf_shock_r <- IRF_conf_shock_r
  res$conf_lvls <- conf_lvls
  res$country <- country
  res$month <- month
  res$year <- year
  res$periods_start <- periods_start
  return(res)
}