eq_generator <- function(y,x,factors=NULL,countrydums=FALSE,cnames=NULL,oprob = FALSE){
  # DESCRIPTION:
  # creates a formula used in the estimation
  #-------------------------------------------------------------------------------
  #USAGE
  # eq_generator(y,x,factors=NULL,countrydums=0,cnames=NULL,oprob = FALSE)
  #-------------------------------------------------------------------------------
  #INPUT
  # y               explained variable
  # x               main explanatory variables (all the differences, and yield levels)
  # factors         vector of factor numbers (7:24)
  # countrydums     should country fixed effects be included?
  # cnames          vector of country names
  # oprob           boolean if this is an oprob equation (without constant)
  #-------------------------------------------------------------------------------
  #OUTPUT
  # formula string
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
idc = "country"
  
eq <- paste(y," ~ ", x[1],sep="")
if (length(x)>1){
  eq <- paste(c(eq,x[2:length(x)]),sep = " + ",collapse=" + ")
}
if (!is.null(factors)){
  for (i in 1:length(factors)){
    eq <- paste(eq," + facratings",factors[i],sep="")
  }
}
if (countrydums){
  print(paste("Check if United States are correctly written: ",which(cnames == "countryUnited.States")))
  for (i in 1:length(cnames)){
    if (nchar(cnames[i])>nchar(idc)){
      if (grepl(idc,cnames[i])){
        if (!grepl("countryUnited.States",cnames[i])){eq <- paste(eq," + ",cnames[i],sep="")}
      }
    }
  }
}

#deduct constant in case of linear model as it conflicts with the penalty matrix
if (!oprob){eq <- paste(eq," - 1", sep="")}
return(eq)
}