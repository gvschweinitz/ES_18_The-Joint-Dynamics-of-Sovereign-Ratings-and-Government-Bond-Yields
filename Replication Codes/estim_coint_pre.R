estim_coint_pre <- function(data,fspecs=fspecs,constlong = FALSE){
  # DESCRIPTION:
  # Performs pre-estimation for a linear model or an ordered probit
  #-------------------------------------------------------------------------------
  #USAGE
  # estim_smooth_pre(formula,data,oprob = FALSE, lambda = 1,polr_start = NULL)
  #-------------------------------------------------------------------------------
  #INPUT
  # data            data frame with the naming convention from data_prep
  # fspecs          Structure controlling estimation setup
  # constlong       boolean to indicate if the constant should be included in the long-run relation
  #-------------------------------------------------------------------------------
  #OUTPUT
  #   pre$LHSlin    LHS of linear model (column names)
  #   pre$LHSoprob  LHS of oprob model (column names)
  #   pre$RHSlin    RHS of linear model (column names)
  #   pre$RHSoprob  RHS of oprob model (column names)
  #   pre$RHSEC     Error correction part (column names)
  #   pre$Ylin      LHS of linear model (data matrix)
  #   pre$Yoprob    LHS of oprob model (data matrix)
  #   pre$Xlin      RHS of linear model (data matrix)
  #   pre$Xoprob    RHS of oprob model (data matrix)
  #   pre$XEC       Error correction part (data matrix)
  #   pre$coeffs_lin    first coefficients linear model
  #   pre$coeffs_oprob  first coefficients oprob model
  #   pre$coeffs_EC     first coefficients error correction part
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz  
  
  library(plyr)
  library(MASS)
  
  contemp <- fspecs$contemp
  cdums <- fspecs$cdums
  lags_dyields_y <- fspecs$ylags_y
  lags_dyields_r <- fspecs$ylags_r
  lags_dratings_y <- fspecs$rlags_y
  lags_dratings_r <- fspecs$rlags_r
  factors <- fspecs$factors
  rasym <- fspecs$asym
  yasym <- fspecs$asym

  data<-cuttolags(data,max(lags_dyields_y,lags_dyields_r),max(lags_dratings_y,lags_dratings_r))
  print(dim(data))
  if (cdums){
    # country fixed effects to be added to eccols (included in longrun-relation) or to lagterm (shortrun relation)
    cnames <- colnames(data)
    cpos <- which(grepl("country",cnames) & (nchar(cnames)>7))
    cpos <- cpos[-length(cpos)]       #no fixed effect for United States
    cadd <- cnames[cpos]
  }
  else{
    cadd <- NULL
  }
  
  formula_lin <- eq_generator_bylag(data,fspecs=fspecs,oprob=FALSE, use_EC = FALSE, setconst=TRUE,long=TRUE)
  formula_oprob <- eq_generator_bylag(data,fspecs=fspecs,oprob=TRUE, use_EC = FALSE, setconst=TRUE,long=TRUE)
#   formula_lin <- eq_generator_bylag(data,lags_dyields_y,lags_dratings_y,oprob = FALSE, rasym = rasym, yasym = yasym,factors = factors,cdums = cdums, use_EC = FALSE, setconst=TRUE,long=TRUE)
#   formula_oprob <- eq_generator_bylag(data,lags_dyields_r,lags_dratings_r,oprob = TRUE, rasym = rasym, yasym = yasym,factors = factors,cdums = cdums, use_EC = FALSE, setconst=TRUE,long=TRUE)
  eq_lin <- as.formula(formula_lin$formula)
  eq_oprob <- as.formula(formula_oprob$formula)
  
  if (constlong){
    # constant (and country fixed effects) estimated in long-run relationship
    eccols <- c("lyields","const",paste("facratings",factors,sep=""),cadd)
    eccols_oprob <- eccols[-pmatch("const",eccols)]
    exog_add <- NULL
  }
  else{
    # constant (and country fixed effects) estimated in short-run relationship
    eccols <- c("lyields",paste("facratings",factors,sep=""))
    eccols_oprob <- eccols
    exog_add <- cadd
  }
  NEC <- length(eccols)
  exog_y <- c(formula_lin$lagterm,exog_add)
  exog_r <- c(formula_oprob$lagterm,exog_add)
  
  res_lin <- lm(eq_lin,data,model=TRUE)
  coeffs_lin <- res_lin$coefficients[exog_y]
#   if (max(diff(coeffs_EC))>40){coeffs_EC <- rep(mean(coeffs_EC),length(coeffs_EC))}
  
  res_oprob <- polr(eq_oprob,data,method="probit")
  if (constlong) {
    coeffs_oprob <- c(res_oprob$coefficients[exog_r],res_oprob$zeta+res_lin$coefficients["const"])
    # use only one threshold as they need to be equally distributed around the constant
    coeffs_oprob <- coeffs_oprob[-(length(coeffs_oprob)-1)]
  }
  else{
    coeffs_oprob <- c(res_oprob$coefficients[exog_r],res_oprob$zeta)
  }
  coeffs_EC_lin <- res_lin$coefficients[eccols]
  coeffs_EC_lin <- coeffs_EC_lin[2:length(coeffs_EC_lin)]/(-coeffs_EC_lin[1])
  coeffs_EC_oprob <- res_oprob$coefficients[eccols_oprob]
  coeffs_EC_oprob <- coeffs_EC_oprob[2:length(coeffs_EC_oprob)]/(-coeffs_EC_oprob[1])
  print(coeffs_EC_lin)
  print(coeffs_EC_oprob)
  coeffs_EC <- (coeffs_EC_lin + coeffs_EC_oprob)/2
  rows <- which(rowSums(is.na(data[,c(as.character(eq_lin)[2],as.character(eq_oprob)[2],exog_y,exog_r,eccols)]))==0)  

  pre <- NULL
  pre$LHSlin <- as.character(eq_lin)[2]
  pre$LHSoprob <- as.character(eq_oprob)[2]
  pre$RHSlin <- exog_y
  pre$RHSoprob <- exog_r
  pre$RHSEC <- eccols
  pre$Ylin <- as.matrix(data[rows,pre$LHSlin])
  pre$Yoprob <- as.matrix(data[rows,pre$LHSoprob])
  pre$Xlin <- as.matrix(data[rows,pre$RHSlin[2:length(pre$RHSlin)]])
  pre$Xoprob <- as.matrix(data[rows,pre$RHSoprob[2:length(pre$RHSoprob)]])
  pre$XEC <- as.matrix(data[rows,pre$RHSEC])
  
  pre$coeffs_lin <- coeffs_lin
  pre$coeffs_oprob <- coeffs_oprob
  pre$coeffs_EC <- coeffs_EC
  print(coeffs_lin)
  print(coeffs_oprob)
  print(coeffs_EC)
  return(pre)
}