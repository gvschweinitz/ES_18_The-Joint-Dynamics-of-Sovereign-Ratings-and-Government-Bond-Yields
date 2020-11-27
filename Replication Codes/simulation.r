simulation <- function(coefficients_y, coefficients_r,fspecs=fspecs, start = c(rating = 24, yield = 6), y_shock = 0, empirical_changes = NULL, periods = 100, output = NULL, details = FALSE, countries = NA) {
  #DESCRIPTION:
  # simulates data using estimated coefficients and shocks. May also plot simulated data
  #-------------------------------------------------------------------------------
  #USAGE
  # simulation(coefficients_y, coefficients_r, powers = 1, rlags_r = 1, ylags_r = 1, rlags_y = 1, ylags_y =  1, asym = FALSE, start = c(rating = 24, yield = 6), y_shock = 0, empirical_changes = NULL, yield_first = TRUE, periods = 100, output = NULL, details = FALSE)
  #-------------------------------------------------------------------------------
  #ARGUMENTS
  # coefficients_y  coefficients of the yield equation
  # coefficients_r  coefficients of the rating equation
  # fspecs          Structure controlling estimation setup
  # start           starting values or scenario (from extract_scenBS)
  # y_shock         Time series of shocks to yield equation (if 0, drawn randomly from normal distribution)
  # empirical_changes matrix giving the empirical distribution of rating changes in the data
  # periods         number of periods to be simulated
  # output          one of "plot", "table" or "last"
  # details         boolean if long results should be returned
  # countries       if not NA, add country dummies
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   results$table     table with simulated ratings and yields
  #   results$X         X matrix for yield equation
  #   results$X_ordinal X matrix for oprob equation
  #   results$Y         Y vector for yield equation
  #   results$Y_ordinal Y vector for oprob equation
  #   results$r_change  vector of rating changes
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  ylags_y <- fspecs$ylags_y
  ylags_r <- fspecs$ylags_r
  rlags_y <- fspecs$rlags_y
  rlags_r <- fspecs$rlags_r
  powers <- fspecs$p_yields
  asym <- fspecs$asym
  cdums <- fspecs$cdums
  final <- fspecs$factors
  yield_first <- fspecs$yield_first
  
  if (!is.list(start)) {scenario = FALSE} else {scenario = TRUE}
  if (length(y_shock)==1) {
    if (y_shock>0) {y_shock <- rnorm(periods,0,y_shock)}else {y_shock <- mat.or.vec(periods,1)}
  } else {
    periods = length(y_shock)
  }
  if(scenario) {
    start_p <- dim(start$yield_X)[1]
    periods = start_p + periods
    y_shock=c(c(mat.or.vec(start_p,1)),c(y_shock))
  } 
  
  #Set shocks
  r_shock <- rnorm(periods,0,1)
  
  #Prepare matrices for yield equation
  X <- mat.or.vec(periods,length(coefficients_y))
  colnames(X) <- names(coefficients_y)
  X[,1+powers] <- 1
  Y <- mat.or.vec(periods,1)
  
  #Prepare matrices for ratings equation
  N_op_coeffs <- length(coefficients_r)
  X_ordinal <- mat.or.vec(periods,N_op_coeffs-2)
  colnames(X_ordinal) <-  names(coefficients_r)[1:(N_op_coeffs-2)]
  Y_ordinal <- mat.or.vec(periods,1)
  r_change <- mat.or.vec(periods,1)
  thresholds <- coefficients_r[(N_op_coeffs-1):N_op_coeffs]
  rating_coeffs <- coefficients_r[1:(N_op_coeffs-2)]  
  
  #Preparing matrices including country dummies
  if (!is.na(countries)) {
    for (c in 1:length(countries)){
      X[,grepl(countries[c],colnames(X))]<-1  
      X_ordinal[,grepl(countries[c],colnames(X_ordinal))]<-1
    }
  }
  #Prepare results table
  table <- mat.or.vec(periods,3)
  colnames(table) <- c('Period','Yield','Rating')
  
  if (!scenario) {
    table[1,1] <- 1
    table[1,2] <- start["yield"]
    table[1,3] <- start["rating"]
    start_period = 2
  } else {
    start_period = start_p +1 
    table[1:start_p,"Period"] = 1:start_p
    table[1:start_p,"Yield"] = start$yield
    table[1:start_p,"Rating"] = start$ratings
    X[1:start_p,] <- as.matrix(start$yield_X)
    X_ordinal[1:start_p,] <- start$rating_X
    r_change[1:start_p] <- start$r_change
    Y[1:start_p] <- start$y_change
  }

  base <- substr(names(coefficients_y),1,10)
  smooth_coeffs_y <- which(base=="facratings")
  numbers_yields <- as.numeric(substr(names(coefficients_y)[smooth_coeffs_y],11,100))
  base <- substr(names(coefficients_r),1,10)
  smooth_coeffs_r <- which(base=="facratings")
  numbers_ratings <- as.numeric(substr(names(coefficients_r)[smooth_coeffs_r],11,100))
 
  yly <- length(ylags_y)
  ylr <- length(ylags_r)
  rly <- length(rlags_y)
  rlr <- length(rlags_r)
  
  rating_list <- function(){
    #### RATINGS ####
    
    # Set rating dummies
    active_dummy <- smooth_coeffs_r[numbers_ratings<=table[t-1,"Rating"]]
    X_ordinal[t,active_dummy] <- 1
    
    # Set rating
    for (p in 1:powers) {X_ordinal[t,p] <- table[t-1,"Yield"]^p}  
    if (!asym) {
      # Set yield differences
      for (l in 1:rly) {if(t-rlags_y[l]>0){X_ordinal[t,powers+l] = Y[t-rlags_y[l]]}}
      # Set rating differences
      for (l in 1:rlr) {if(t-rlags_r[l]>0){X_ordinal[t,powers+rly+l] = r_change[t-rlags_r[l]]}}
    } else {
      # Set yield differences
      for (l in 1:rly) {if(t-rlags_y[l]>0){
        if(Y[t-rlags_y[l]]>0){
          X_ordinal[t,powers+l] = Y[t-rlags_y[l]]
        }}}
      for (l in 1:rly) {if(t-rlags_y[l]>0){if(Y[t-rlags_y[l]]<0){
        X_ordinal[t,powers+rly+l] = Y[t-rlags_y[l]]
      }}}
      # Set rating differences
      for (l in 1:rlr) {if(t-rlags_r[l]>0){if(r_change[t-rlags_r[l]]>0){X_ordinal[t,powers+2*rly+l] = r_change[t-rlags_r[l]]}}}
      for (l in 1:rlr) {if(t-rlags_r[l]>0){if(r_change[t-rlags_r[l]]<0){X_ordinal[t,powers+2*rly+rlr+l] = r_change[t-rlags_r[l]]}}}
    }
    
    # Compute Y(t) 
    Y_ordinal[t] <- sum(X_ordinal[t,]*rating_coeffs)+r_shock[t]
    if (Y_ordinal[t]<thresholds[1]) {r_change[t]=-1} else {if(Y_ordinal[t]>thresholds[2]){r_change[t]=+1}}
    
    if ((is.null(empirical_changes)) | (r_change[t]==0)) {
      current_rating <- table[t-1,"Rating"]+r_change[t]
    } else {
      if (r_change[t]>0) {
        current_rating <- table[t-1,"Rating"]+empirical_changes[which(runif(1)<empirical_changes[,2])[1],1]
      } else {
        current_rating <- table[t-1,"Rating"]-empirical_changes[which(runif(1)<empirical_changes[,3])[1],1]
      }
    }
    if (current_rating<6) {
      current_rating=6
      r_change[t] = 0
    }
    if (current_rating>24) {
      current_rating=24
      r_change[t]=0
    }
    
    rating_list <- NULL
    rating_list$current_rating <- current_rating
    rating_list$r_change <- r_change
    rating_list$Y_ordinal <- Y_ordinal
    rating_list$X_ordinal <- X_ordinal
    return(rating_list)
  }
  yield_list <- function(){
    #### YIELDS ####
    
    # Set rating dummies
    active_dummy <- smooth_coeffs_y[numbers_yields<=table[t-1,"Rating"]]
    X[t,active_dummy] <- 1
    
    # Set interest level
    for (p in 1:powers) {X[t,p] <- table[t-1,"Yield"]^p}
    
    if (!asym) {
      # Set yield differences
      for (l in 1:yly) {if(t-ylags_y[l]>0){X[t,1+powers+l] = Y[t-ylags_y[l]]}}
      # Set rating differences
      for (l in 1:ylr) {if(t-ylags_r[l]>0){X[t,1+powers+yly+l] = r_change[t-ylags_r[l]]}}
    } else {
      # Set yield differences
      for (l in 1:yly) {if(t-ylags_y[l]>0){if(Y[t-ylags_y[l]]>0){
        X[t,1+powers+l] = Y[t-ylags_y[l]]}}}
      for (l in 1:yly) {if(t-ylags_y[l]>0){if(Y[t-ylags_y[l]]<0){
        X[t,1+powers+yly+l] = Y[t-ylags_y[l]]}}}
      # Set rating differences
      for (l in 1:ylr) {if(t-ylags_r[l]>0){if(r_change[t-ylags_r[l]]>0){X[t,1+powers+2*yly+l] = r_change[t-ylags_r[l]]}}}
      for (l in 1:ylr) {if(t-ylags_r[l]>0){if(r_change[t-ylags_r[l]]<0){X[t,1+powers+2*yly+ylr+l] = r_change[t-ylags_r[l]]}}}
    }
    
    # Compute Y(t)
    Y[t] <- sum(X[t,]*coefficients_y) + y_shock[t]
    current_yield <- table[t-1,"Yield"]+Y[t]
    
    yield_list <- NULL
    yield_list$current_yield <- current_yield
    yield_list$y_shock <- y_shock
    yield_list$Y <- Y
    yield_list$X <- X
    return(yield_list)
  }

  for (t in start_period:periods) {

    table[t,1] = t
    
    if (yield_first){
      yl <- yield_list()
      table[t,"Yield"] <- yl$current_yield
      y_shock <- yl$y_shock
      Y <- yl$Y
      X <- yl$X
      
      rl <- rating_list()
      table[t,"Rating"] <- rl$current_rating
      r_change <- rl$r_change
      Y_ordinal <- rl$Y_ordinal
      X_ordinal <- rl$X_ordinal
    }
    else{
      rl <- rating_list()
      table[t,"Rating"] <- rl$current_rating
      r_change <- rl$r_change
      Y_ordinal <- rl$Y_ordinal
      X_ordinal <- rl$X_ordinal
      
      yl <- yield_list()
      table[t,"Yield"] <- yl$current_yield
      y_shock <- yl$y_shock
      Y <- yl$Y
      X <- yl$X
    }
  }
  if (!is.null(output)) {
    switch(output,
           plot = {plot(table[,'Rating'])
                   if(scenario){lines(c(start_p,start_p),c(min(table[,"Rating"]),max(table[,"Rating"])))}
                   if(scenario){lines(c(0,periods),c(table[start_p,"Rating"],table[start_p,"Rating"]))}},
           table = {print(table)})
           last = {print(table[periods,])}
  }

  if (details == FALSE) {
    return(table)
  } else {
    results <- NULL
    results$table <- table
    results$X <- X
    results$X_ordinal <- X_ordinal
    results$Y <- Y
    results$Y_ordinal <- Y_ordinal
    results$r_change <- r_change
    return(results)
  }
}