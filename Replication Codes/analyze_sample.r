analyze_sample <- function(data,yield_results, crossid = "country",dateid = "dateid",maxlag = 2) {
  # DESCRIPTION:
  # prepares information on observations for wild bootstrap
  #-------------------------------------------------------------------------------
  #USAGE
  # analyze_sample(data,yield_results, crossid = "country",dateid = "dateid",maxlag = 2)
  #-------------------------------------------------------------------------------
  #INPUT
  # data          data frame with data with the naming convention from data_prep
  # yield_results results from yield estimation
  # crossid       name of country column (N cross sections)
  # dateid        name of date column (T observations)
  # maxlag        maximum number of lags to be accounted for -> additional lags of residuals may have to be deleted
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
  
  
  ids <- unique(data[,crossid])
  dates <- as.character(unique(data[,dateid]))
  
  #works for old and new estimation style
  residuals <- yield_results$res_smooth$residuals
  orig_line <- as.numeric(names(residuals))
  
  N <- length(ids)
  T <- length(dates)
  
  date_val <- mat.or.vec(T,1)
  
  for (t in 1:T) {
    s <- strsplit(dates[t],split = "/")
    month <- as.numeric(s[[1]][1])
    year <- as.numeric(s[[1]][3])
    date_val[t] <- year+month/13
  }
  dates <- dates[order(date_val)]
  
  resid_mat <- mat.or.vec(T,N)+NA
  colnames(resid_mat) <- ids
  rownames(resid_mat) <- dates
  loc_mat <- resid_mat
  
  
  for (id in ids) {
    for (d in dates) {
      i <- which((data[,crossid]==id) & (data[,dateid]==d))
      i_old <- i
      i <- which(orig_line==i)
      if (length(i)>0) {
        resid_mat[d,id] <- residuals[i]
        loc_mat[d,id] <- i_old
      }
    }
  }
  
  loc_copy <- loc_mat
  loc_original <- loc_mat
  loc_copy[!is.na(loc_copy)] <- 0
  loc_mat[1:maxlag,]=NA
  loc_mat[(maxlag+1):T,] <- loc_mat[(maxlag+1):T,]+loc_copy[1:(T-maxlag),]
  loc_vec <- as.matrix(c(loc_mat))
  loc_vec <- loc_vec[-which(is.na(loc_vec))]
  
  results <- NULL
  results$resid_mat <- resid_mat
  results$loc_original <- loc_original # all ids
  results$loc_mat <- loc_mat
  results$loc_vec <- loc_vec  # eligible "seeds"
  results$T <- T
  results$N <- N
  return(results)
}