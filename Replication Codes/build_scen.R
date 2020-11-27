build_scen <- function(yields,ratings,coefficients_y,coefficients_r,fspecs=fspecs){
  #DESCRIPTION:
  # constructs scenario (to be used in "simulation.R") from yield and rating sequence
  #-------------------------------------------------------------------------------
  #USAGE
  # build_scen(yields,ratings,coefficients_y,coefficients_r,fspecs=fspecs)
  #-------------------------------------------------------------------------------
  #INPUT
  # yields            Initial yield levels to construct scenario
  # ratings           Initial rating levels to construct scenario
  # coefficients_y    coefficient vector of yield equation. Only needed to know variables to be constructed
  # coefficients_r    coefficient vector of rating equation. Only needed to know variables to be constructed
  # fspecs            Structure controlling estimation setup
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   yield           yields
  #   ychange         first difference of yields (LHS variable in yield equation)
  #   ratings         ratings
  #   rchange         upgrade or downgrade (LHS variable in rating equation)
  #   yield_X         RHS variables for yield equation
  #   rating_X        RHS variables for rating equation
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  ylags_y <- fspecs$ylags_y
  ylags_r <- fspecs$ylags_r
  rlags_y <- fspecs$rlags_y
  rlags_r <- fspecs$rlags_r
  powers <- fspecs$p_yields
  asym <- fspecs$asym
  
  if (length(yields)==1){
    return(c(rating=ratings,yield=yields))
  }
  
  y_change <- c(0,diff(yields))
  r_change <- c(0,sign(diff(ratings)))
  T <- length(yields)
  yly <- length(ylags_y)
  ylr <- length(ylags_r)
  rly <- length(rlags_y)
  rlr <- length(rlags_r)
  
  base <- substr(names(coefficients_y),1,10)
  smooth_coeffs_y <- which(base=="facratings")
  numbers_yields <- as.numeric(substr(names(coefficients_y)[smooth_coeffs_y],11,100))
  base <- substr(names(coefficients_r),1,10)
  smooth_coeffs_r <- which(base=="facratings")
  numbers_ratings <- as.numeric(substr(names(coefficients_r)[smooth_coeffs_r],11,100))
  
  X <- mat.or.vec(T,length(coefficients_y))
  X_ordinal <- mat.or.vec(T,length(coefficients_r)-2)
  colnames(X) <- names(coefficients_y)
  colnames(X_ordinal) <- names(coefficients_r)[1:(length(coefficients_r)-2)]
  
  for (p in 1:powers) {
    X[2:T,p] <- yields[1:(T-1)]^p
    X_ordinal[2:T,p] <- yields[1:(T-1)]^p
  }
  X[,1+powers] <- 1
  
  for (t in 2:T){
    active_dummy <- smooth_coeffs_y[numbers_yields<=ratings[t-1]]
    X[t,active_dummy] <- 1
    active_dummy <- smooth_coeffs_r[numbers_ratings<=ratings[t-1]]
    X_ordinal[t,active_dummy] <- 1
  }
  
  if (!asym) {
    # Set yield differences
    for (l in 1:yly) {
      if(T-l>0){
        X[(ylags_y[l]+1):T,1+powers+l] = y_change[1:(T-ylags_y[l])]
      }
    }
    for (l in 1:ylr) {
      if(T-l>0){
        X_ordinal[(ylags_r[l]+1):T,powers+l] = y_change[1:(T-ylags_r[l])]
      }
    }
    # Set rating differences
    for (l in 1:rly) {
      if(T-l>0){
        X[(rlags_y[l]+1):T,1+powers+yly+l] = r_change[1:(T-rlags_y[l])]
      }
    }
    for (l in 1:rlr) {
      if(T-l>0){
        X_ordinal[(rlags_r[l]+1):T,powers+ylr+l] = r_change[1:(T-rlags_r[l])]
      }
    }
  } 
  else {
    # Set yield differences
    for (l in 1:yly) {
      if(T-l>0){
        X[(ylags_y[l]+1):T,1+powers+l] = pmax(y_change[1:(T-ylags_y[l])],0)
        X[(ylags_y[l]+1):T,1+powers+yly+l] = pmin(y_change[1:(T-ylags_y[l])],0)
      }
    }
    for (l in 1:ylr) {
      if(T-l>0){
        X_ordinal[(ylags_r[l]+1):T,powers+l] = pmax(y_change[1:(T-ylags_r[l])],0)
        X_ordinal[(ylags_r[l]+1):T,powers+ylr+l] = pmin(y_change[1:(T-ylags_r[l])],0)
      }
    }
    # Set rating differences
    for (l in 1:rly) {
      if(T-l>0){
        X[(rlags_y[l]+1):T,1+powers+2*yly+l] = pmax(r_change[1:(T-rlags_y[l])],0)
        X[(rlags_y[l]+1):T,1+powers+2*yly+rly+l] = pmin(r_change[1:(T-rlags_y[l])],0)
      }
    }
    for (l in 1:rlr) {
      if(T-l>0){
        X_ordinal[(rlags_r[l]+1):T,powers+2*ylr+l] = pmax(r_change[1:(T-rlags_r[l])],0)
        X_ordinal[(rlags_r[l]+1):T,powers+2*ylr+rlr+l] = pmin(r_change[1:(T-rlags_r[l])],0)
      }
    }
  }
  
  res <- NULL
  res$yield <- yields
  res$y_change <- c(0,diff(yields))
  res$ratings <- ratings
  res$r_change <- c(0,sign(diff(ratings)))
  res$yield_X <- X
  res$rating_X <- X_ordinal
  
  return(res)
}