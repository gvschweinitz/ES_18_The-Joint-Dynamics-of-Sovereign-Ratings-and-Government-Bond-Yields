estim_smooth_pre <- function(formula,data,oprob = FALSE,lambda=1,polr_start = NULL){
  # DESCRIPTION:
  # Performs pre-estimation for a linear model or an ordered probit
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_smooth_pre(formula,data,oprob = FALSE, lambda = 1,polr_start = NULL)
  #-------------------------------------------------------------------------------
  #ARGUMENTS
  # formula         formula to be estimated
  # data            data frame with the naming convention from data_prep
  # oprob           boolean to indicate ordered probit
  # lambda          smoothing parameter for facratings (only needed to check if linear or full estimation is to be done)
  # polr_start      starting values for the polr estimation (lin or general, in turn used as starting values for LL)
  #                   values take the form c(coefficients,zeta)
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   pre$X         matrix with X values for later estimation (cut to availability)
  #   pre$Y         matrix with Y values for later estimation (cut to availability). For oprob, values are 1,2,3
  #   pre$RHS       vector of strings denoting RHS variables
  #   pre$LHS       string denoting LHS variable
  #   pre$coeffs_start  vector of starting coefficients
  #   pre$coeffs_lin    vector of coefficients from oprob with linearity restriction (for comparison)
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  library(plyr)
  library(MASS)
  
  lambdathresh <- 1000^3
  if (!is.formula(formula)) {
    formula <- as.formula(formula)
  }
  
  LHS <- as.character(formula)[2]
  RHS <- strsplit(as.character(formula)[3],' - ',fixed = TRUE)[[1]][1]     #delete -1
  RHS <- strsplit(RHS,' + ',fixed = TRUE)[[1]]
  cut_coeffs <- grepl("\n",RHS)
  RHS[cut_coeffs] <- substr(RHS[cut_coeffs],6,nchar(RHS[cut_coeffs]))

  if (oprob) {
    #oprob with two cases: linear estimation in case of high lambda, nonlinear otherwise
    if (lambda<lambdathresh){
      res_lin <- NULL
      if (is.null(polr_start)){res_start <- polr(formula, data, method = "probit", model = TRUE)}
      else {res_start <- polr(formula, data, start=polr_start, method = "probit", model = TRUE)}
      rows <- rownames(res_start$model)
      coeffs_lin <- NULL
      coeffs_start <- c(res_start$coefficients,res_start$zeta)
    }
    else{
      base <- substr(RHS,1,10)
      RHS_lin <- c(RHS[which(base!="facratings")],"lgratings2")
      ratpos <- which(RHS_lin=="lgratings2")
      if (sum(grepl("lgratings2",colnames(data)))==0){data$lgratings2 <- data$lgratings-6}
      formula_lin <- as.formula(paste(LHS,"~",paste(RHS_lin,collapse=" + ")))
      if (is.null(polr_start)){res_lin <-polr(formula_lin,data,method = "probit",model=TRUE)}
      else {
        
        if (length(polr_start) == length(RHS)+2){
          polr_start <- c(polr_start[which(base!="facratings")],
                          mean(polr_start[which(base=="facratings")]),
                          polr_start[(length(RHS)+1):(length(RHS)+2)])
        }
        res_lin <- polr(formula_lin, data, start=polr_start, method = "probit", model = TRUE)
      }
      rows <- rownames(res_lin$model)
      coeffs_lin <- c(res_lin$coefficients,res_lin$zeta)
      
      #expand linear coefficients to nonlinear case
      coeffs_start <- rep(0,length(RHS))
      names(coeffs_start) <- RHS
      coeffs_start[RHS[which(base!="facratings")]] <- coeffs_lin[RHS_lin[which(RHS_lin!="lgratings2")]]
      coeffs_start[RHS[which(base=="facratings")]] <- coeffs_lin[RHS_lin[which(RHS_lin=="lgratings2")]]
      coeffs_start <- c(coeffs_start,res_lin$zeta)
    }
  }
  else{
    res_start <- lm(formula, data, x = TRUE)
    rows <- rownames(res_start$x)
    res_start$x <- NULL
    coeffs_start <- res_start$coefficients
    coeffs_lin <- NULL
  }

  pre <- NULL
  pre$X <- as.matrix(data[rows,RHS])
  pre$Y <- as.matrix(data[rows,LHS])
  if(oprob){pre$Y <- as.numeric(pre$Y)+2}     #changed from -1,0,1 to 1,2,3
  pre$RHS <- RHS
  pre$LHS <- LHS
  pre$coeffs_start <- coeffs_start
  pre$coeffs_lin <- coeffs_lin
  return(pre)
}