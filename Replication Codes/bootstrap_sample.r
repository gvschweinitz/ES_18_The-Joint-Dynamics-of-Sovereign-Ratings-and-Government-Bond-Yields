bootstrap_sample <- function(sample_info, empirical_changes,data,res_y,res_r,fspecs=fspecs) {
  #DESCRIPTION:
  # Creates bootstrapped data set using a wild bootstrap.
  # For every country, a starting scenario (maxlag sequence of yields and ratings) is drawn at random
  # Starting from this scenario, the development of this country is simulated forward using the estimated coefficients and a wild bootstrap of original estimation errors
  # The wild bootstrap accounts for cross-country correlation (same bootstrap factor applied to all countries at a given time)
  # WARNING: THE SCENARIOS DO NOT ACCOUNT FOR DIFFERENT COUNTRIES AT THE MOMENT. IMPLICATION: bootstrap of country dummies is misleading
  #-------------------------------------------------------------------------------
  #USAGE
  # bootstrap_sample(sample_info, empirical_changes,data,res_y,res_r,fspecs=fspecs)
  #-------------------------------------------------------------------------------
  #INPUT
  # sample_info       Result from analyze_sample.R
  # empirical_changes Result from extract_changes.R
  # data              data frame with data with the naming convention from data_prep
  # res_y             results from yield estimation
  # res_r             results from rating estimation
  # fspecs            Structure controlling estimation setup
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   resid_mat     TxN matrix of estimation residuals
  #   log_original  TxN matrix of row indices wrt original data
  #   loc_mat       TxN matrix of row indices wrt original data after accounting for lags
  #   loc_vec       vector with all non-NaN-elements from loc_mat
  #   T             # of different crossid
  #   N             # of different dateid
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  ylags_y = fspecs$ylags_y
  ylags_r = fspecs$ylags_r
  rlags_y = fspecs$rlags_y
  rlags_r = fspecs$rlags_r
  asym = fspecs$asym
  cdums <- fspecs$cdums
  rat_var = fspecs$idratings
  yield_var = fspecs$idyields_final  
  
  pre_y <- res_y$pre
  pre_r <- res_r$pre

  y_dep_name <- "dyields0"
  r_dep_name <- "dratings0_fac"
  
  y_shocks <- sample_info$resid_mat
  T <- length(y_shocks[!is.na(y_shocks)])

  block <- mat.or.vec(T,length(res_y$res_smooth$coefficients)+1)
  block_o <- mat.or.vec(T,length(res_r$res_smooth$coefficients)-2+1)
  colnames(block) <- c(y_dep_name,names(res_y$res_smooth$coefficients))
  colnames(block_o) <- c("RDEP",names(res_r$res_smooth$coefficients)[1:(length(res_r$res_smooth$coefficients)-2)])
  dist <- c(-sqrt(1.5),-1,-sqrt(0.5),sqrt(0.5),1,sqrt(1.5))
  N_dist = length(dist)
  mult <- dist[sample(1:N_dist,sample_info$T,replace=TRUE)]
  
  for (t in 1:sample_info$T) {
   y_shocks[t,] <- y_shocks[t,]*mult[t]
  }
  
  max_lag <- max(rlags_r,ylags_r,rlags_y,ylags_y)
  tt <-0
  for (c in 1:sample_info$N) {
    #print(c)
    if (cdums){country <- colnames(y_shocks)[c]} else{country<-NA}
    
    # Scenario is not country-specific at the moment
    check <- FALSE
    while (!check){
      start <- sample_info$loc_vec[sample(1:length(sample_info$loc_vec),1)]

      scen <- extract_scenBS((start-max_lag+1):start,data,pre_y,pre_r,rat_var = rat_var, yield_var = yield_var, country = country)
      check <- (sum(is.na(scen$yield))==0) & 
                (sum(is.na(scen$ratings))==0) &
                (sum(is.na(scen$yield_X))==0) &
                (sum(is.na(scen$rating_X))==0) &
                (sum(is.na(scen$r_change))==0) &
                (sum(is.na(scen$y_change))==0)
      if (!check){
        print("draw again")
      }
    }
    y_shock <- y_shocks[,c]
    y_shock <- y_shock[-which(is.na(y_shock))]
    t <- length(y_shock)

    sim_c <- simulation(res_y$res_smooth$coefficients, res_r$res_smooth$coefficients,fspecs=fspecs, start = scen, y_shock = y_shock, empirical_changes = empirical_changes, output = NULL, details = TRUE, countries = country)
    cp_block <- cbind(sim_c$Y[(max_lag+1):(t+max_lag)],sim_c$X[(max_lag+1):(t+max_lag),])
    cp_block_o <- cbind(sim_c$r_change[(max_lag+1):(t+max_lag)],sim_c$X_ordinal[(max_lag+1):(t+max_lag),])
    block[(tt+1):(tt+t),] <- cp_block
    block_o[(tt+1):(tt+t),] <- cp_block_o
    tt <- tt+t
  }
  
  block <- as.data.frame(block)
  block_o <- as.data.frame(block_o)
  block_o$dratings0_fac <- as.factor(block_o$RDEP)
  facratings_i = substr(colnames(block_o),1,10)=="facratings"
  block_o$lgratings2 <- rowSums(block_o[,facratings_i])

  result <- NULL
  result$data_y <- block
  result$data_r <- block_o 

  return(result)
}