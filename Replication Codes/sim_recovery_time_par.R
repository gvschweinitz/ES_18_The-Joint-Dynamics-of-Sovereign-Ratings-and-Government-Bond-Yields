sim_recovery_time_par <- function(res_yields, res_ratings,fspecs=fspecs, start=c(rating = 24, yield = 6), empirical_changes = NULL, periods = 100, determ = FALSE) {
  #DESCRIPTION:
  # Calculates the first occurrence of every rating in months on the recovery path from a low rating/yield combination. Returns results from 5000 simulations
  #-------------------------------------------------------------------------------
  #USAGE
  # sim_recovery_time_par(res_yields, res_ratings,fspecs=fspecs, start=c(rating = 24, yield = 6), empirical_changes = NULL, periods = 100, determ = FALSE)
  #-------------------------------------------------------------------------------
  #INPUT
  # res_yields        results from yield estimation
  # res_ratings       results from rating estimation
  # fspecs            Structure controlling estimation setup
  # start             Initial rating and yield
  # empirical_changes Result from extract_changes.R
  # periods           Simulation horizon.
  # determ            Boolean. FALSE: usual yield shocks may occur over the IRF duration. TRUE: upgrade/downgrade probability is the only source of uncertainty in simulation
  #-------------------------------------------------------------------------------
  #OUTPUT
  # time              5000 x #facratings matrix, containing month of first occurrence of a given rating in a given simulation draw
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  library(foreach)
  library(doParallel)
  print("opening clusters")
  cl<-makeCluster(16, outfile="")
  registerDoParallel(cl)
  functionlist <- c("simulation")
  
  cdums <- fspecs$cdums
  final <- fspecs$factors
  
  reps = 5000
  initialize = 20
  
  if (cdums) {countries="United.States"}
  else{countries=NA}
  
  coefficients_r <- res_ratings$res_smooth$coefficients
  coefficients_y <- res_yields$res_smooth$coefficients
  if (determ) {y_shock <- 0}
  else {y_shock <- sd(res_yields$res_smooth$residuals[1:length(res_yields$pre$Y)])}
  
  time <- mat.or.vec(reps,length(final))
  colnames(time) <- final
  
  lc = 1
  
  calc_time <- function(i){
    r <- simulation(coefficients_y = coefficients_y, coefficients_r = coefficients_r,fspecs=fspecs, start = start, y_shock = y_shock, empirical_changes = empirical_changes, periods = periods,countries=countries)
    time_temp <- mat.or.vec(1,length(final))
    colnames(time_temp) <- final
    for (f in 1:length(final)) {
      t = which(r[,3]>=final[f])
      if (length(t)>0) {time_temp[f] <- t[1]} else {time_temp[f] = Inf}
    }
    return(time_temp)
  }
  
  time <- foreach (i=1:reps,.combine=rbind,.export=functionlist) %dopar% (calc_time(i))

  return(time)
}