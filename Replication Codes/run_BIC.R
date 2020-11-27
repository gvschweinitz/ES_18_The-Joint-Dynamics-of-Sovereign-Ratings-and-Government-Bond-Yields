run_BIC  <-  function(data,maxlags,fspecs = fspecs,oprob = FALSE,optimcol = 3,optimlambda = TRUE){
  # DESCRIPTION:
  # finds the optimal lag order
  # Data are cut to maximum lag order to make results comparable. Optimality criterion is selected based on optimal lambda
  #-------------------------------------------------------------------------------
  #USAGE
  # run_BIC(data,maxlags,oprob = FALSE,optimcol = 3,cdums = 0,cnames=NULL,optimlambda = TRUE)
  #-------------------------------------------------------------------------------
  #INPUT
  # data            data frame with the naming convention from data_prep
  # maxlags         maximum number of lags in the yield equation
  # oprob           boolean to indicate ordered probit
  # optimcol        column in the result-table after which one should optimize
  # cdums           boolean to indicate if country dummies should be used
  # cnames          names of countries (to identify country dummies)
  # optimlambda     boolean to indicate if lambda should be fully optimized
  #-------------------------------------------------------------------------------
  #OUTPUT
  # out             maxlags x 5 matrix with columns lambda, AICc, BIC, kappa and LL_data
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  rasym <- fspecs$asym
  yasym <- fspecs$asym
  contemp <- fspecs$contemp
  factors <- fspecs$factors
  cdums <- fspecs$cdums
  yield_first <- fspecs$yield_first
  fspecs_run <- fspecs
  
  data <- cuttolags(data,maxlags,maxlags)
  cnames <- colnames(data)
  out <- mat.or.vec(maxlags,5)
  
  for (l in 1:maxlags){
    if (oprob){
      # Run estimation of ordered-probit rating equation
      if (yield_first) {fspecs_run$rlags_y <- c(0:l)}
      else{fspecs_run$rlags_y <- c(1:(l+1))}
      if (!contemp){fspecs_run$rlags_y <- c(1:l)}
      fspecs_run$rlags_r <- c(1:l)
      formula <- eq_generator_bylag(data,fspecs=fspecs_run,oprob=oprob)

      if (optimlambda){
        # optimize lambda?
        res_dyn <- IC_all_dyn(formula=formula, data = data, max_diff = 2^14, oprob = TRUE, func = function(x){2^10*(x-1)},optimcol = optimcol)
        opt_row <- res_dyn[res_dyn[,optimcol]==min(res_dyn[,optimcol]),]
      }
      else {
        # Run regression with full rating smoothing
        res <- estim_smooth(formula=formula,data=data,lambda = 1000^3,oprob = TRUE, calcgr=FALSE, calchess = FALSE, usemulti = FALSE, polr_start = NULL, calcpost=TRUE)
        opt_row <- c(res$lambda,res$post$AICc,res$post$BIC,res$kappa,res$post$LL_out$LL_data)
      }
    }
    else {
      # Run estimation of OLS yield equation, always with lambda optimization
      if (yield_first) {fspecs_run$ylags_r <- c(1:(l+1))}
      else{fspecs_run$ylags_r <- c(0:l)}
      if (!contemp){fspecs_run$ylags_r <- c(1:l)}
      fspecs_run$ylags_y <- c(1:l)
      formula <- eq_generator_bylag(data,fspecs=fspecs_run,oprob=oprob)
      res_dyn <- IC_all_dyn(formula=formula, data = data, max_diff = 0.01, oprob = FALSE, func = function(x){max(2*(x-1),0)},optimcol = optimcol)
      opt_row <- res_dyn[res_dyn[,optimcol]==min(res_dyn[,optimcol]),]
    }
    
    out[l,] <- opt_row
    print(out)
  }
  
  if (optimlambda){colnames(out) <- colnames(res_dyn)}
  else {colnames(out) <- c("lambda","AICc","BIC","kappa","LL_data")}
  
  return(out)
}