extract_scenBS <- function(lines,data,pre_y,pre_r,rat_var = "ratings_f_last", yield_var = "ryields", country = NA) {
  #DESCRIPTION:
  # Extract scenario (starting values) for bootstrap draws
  #-------------------------------------------------------------------------------
  #USAGE
  # extract_scenBS(lines,data,pre_y,pre_r,rat_var = "ratings_f_last", yield_var = "ryields", country = NA)
  #-------------------------------------------------------------------------------
  #INPUT
  # lines         Selection observations of the scenario
  # data          data frame with data with the naming convention from data_prep
  # pre_y         Output of estim_smooth_pre.R for yield estimation
  # pre_r         Output of estim_smooth_pre.R for rating estimation
  # rat_var       name of rating variable
  # yield_var     name of yield variable
  # country       If necessary, add country fixed effect
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   yield       yields of starting scenario
  #   ratings     ratings of starting scenario
  #   yield_X     data of starting scenario for yield estimation
  #   rating_X    data of starting scenario for rating estimation
  #   y_change    yield changes of starting scenario (LHS yield equation)
  #   r_change    rating changes of starting scenario (LHS rating equation)
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  y_select <- data[as.character(lines),c(pre_y$LHS,pre_y$RHS)]
  y_change <- as.matrix(y_select[,1])
  y_select <- as.matrix(y_select[,-1])
  
  r_select <- data[as.character(lines),c(pre_r$LHS,pre_r$RHS)]
  r_select[,pre_r$LHS] <- as.numeric(r_select[,pre_r$LHS])-2
  r_change <- as.matrix(r_select[,1])
  r_select <- as.matrix(r_select[,-1])
  

# Deal with country fixed effects
  if (!is.na(country)){
    y_select[,grepl("country",colnames(y_select))] <- 0
    r_select[,grepl("country",colnames(r_select))] <- 0
    y_select[,grepl(country,colnames(y_select))] <- 1
    r_select[,grepl(country,colnames(r_select))] <- 1
  }
  result <- NULL
  result$yield <- data[lines,yield_var]
  result$ratings <- data[lines,rat_var]
  result$yield_X <- y_select
  result$rating_X <- r_select
  result$r_change <- r_change
  result$y_change <- y_change
  return(result)
}