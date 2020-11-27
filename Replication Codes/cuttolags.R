cuttolags <- function (data, maxlags_y,maxlags_r){
  # DESCRIPTION:
  # creates a consistent dataset for lag-order comparison (putting NA's corresponding to the max. lag-order)
  #-------------------------------------------------------------------------------
  #USAGE
  # cuttolags(data,lags_dyields,lags_dratings)
  #-------------------------------------------------------------------------------
  #INPUT
  # data          data frame with data with the naming convention from data_prep
  # maxlags_y     maximum lag order of differenced yields
  # maxlags_r     maximum lag order of differenced ratings
  #-------------------------------------------------------------------------------
  #OUTPUT
  # data          adjusted data.frame
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  col_y <- paste("dyields",c(0:maxlags_y),sep="")
  col_r <- paste("dratings",c(0:maxlags_r),sep="")
  cntNA <- rowSums(is.na(data[,c(col_y,col_r)]))
  putNA <- c("ryields","lyields","dratings0_fac",col_y,col_r)
  data[cntNA>0,putNA] <- NA
  
  
  return(data)
}