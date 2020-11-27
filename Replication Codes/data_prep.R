data_prep <- function(data_file,lags_dyields,lags_dratings,fspecs){
  #DESCRIPTION
  # data_prep loads data from a csv, calculates real yields and lags of differences, appends country fixed effects, and deletes NA-rows
  #-------------------------------------------------------------------------------
  #USAGE
  # data_prep(data_file,lags_dyields,lags_dratings,fspecs)
  #-------------------------------------------------------------------------------
  #INPUT
  # data_file       name of the file (without .csv)
  # lags_dyields    number of lags of differenced yields
  # lags_dratings   number of lags of differenced ratings
  # fspecs          Structure containing information on
  #   idc             column name of country identifier
  #   idobs           column name of date identifier
  #   idratings       column name of rating information
  #   idyields        column name of yield information
  #   iddeflator      column name of deflator information (in case realyields==TRUE)
  #   asym            generates separate data series for positive and negative values of yields an ratings
  #   p_yields        polynomial degree of yields (nonlinear cointegration relationship - should not be used)
  #   realyields      boolean if real yields should be calculated and used (using the deflator indicated in the file)
  #   real_exp        boolean if (rational) expectations be used to deflate yields: deflate yields at t with deflator at t+12
  #   cdums           boolean if country dummies should be appended
  #   na.rm           boolean if rows with NA in the two key variables should be removed
  #-------------------------------------------------------------------------------
  #OUTPUT
  # data            data.frame (dateid and crossid) with all necessary LHS and RHS variables
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
library(data.table)
library(plyr)
library(dummies)

idc = fspecs$idc
idobs = fspecs$idobs
idratings= fspecs$idratings
idyields=fspecs$idyields
iddeflator = fspecs$iddeflator
asym = fspecs$asym 
p_yields = fspecs$p_yields
realyields = fspecs$realyields
real_exp=fspecs$real_exp
cdums = fspecs$cdums
na.rm = fspecs$na.rm

dta = read.csv(paste(data_file,".csv",sep=""),sep = ",")
dta <- as.data.table(dta)

## PREPARE YIELDS
if (is.factor(dta[[idyields]])) {dta[[idyields]] <- as.numeric(as.character(dta[[idyields]]))}
if (realyields){
  # Calculate real yields from nominal yields and deflator
  if (na.rm){
    data <- dta[!is.na(get(idratings)) & ! is.na(get(idyields)) & !is.na(get(iddeflator)),]
  }else{
    data <- dta
  }
  data[,const:=1]
  
  if (real_exp) {
    data[,ryields := get(idyields)-shift(get(iddeflator),n=12,type="lead"),by=get(idc)]
  }else{
    data[,ryields := get(idyields)-get(iddeflator)]
  }
  
  # Remove real yields at threshold above 200%
  data[ryields>200,ryields := NA]
}else{
  if (na.rm) {
    data <- dta[!is.na(get(idratings)) & ! is.na(get(idyields)),]
  }else{
    data <- dta
  }
  data[,const:=1]
  data[,ryields := get(idyields)]
}

# Add lagged real yields
data[,lyields := shift(ryields,n=1,type="lag"),by=get(idc)]
if (p_yields>1){
  # Add yield polynomial
  for (i in 2:p_yields){
    data[temp := lyields^i]
    setnames(data,"temp",paste("lyields^",i,sep=""))
  }
}

# Add (lags of) yield changes
data[,dyields0 := ryields-shift(ryields,n=1,type="lag"),by=get(idc)]
for (i in 1:lags_dyields){
  data[,paste0("dyields",i) := shift(dyields0,n=i,type="lag"),by=get(idc)]
}

## PREPARE RATINGS
# Create rating groups
data[,gratings := round(get(idratings),digits=0)]
data[,lgratings := shift(gratings,n=1,type="lag"),by=get(idc)]
#adjustment for scaling on 1:18 instead of 7:24
data[,lgratings2 := lgratings-6]
data[lgratings2<=0,lgratings2 := 0]
# Add (cumulative) rating dummies for different rating classes, used for Breitung-Roling smoothing function
data[,facratings := as.factor(lgratings)]
dum <- dummy("facratings",data)
if (length(which(colnames(dum)=="facratingsNA"))>0) {dum <- dum[,-which(colnames(dum)=="facratingsNA")]}
nc <- ncol(dum)
for (i in (nc-1):1) {dum[dum[,i+1]>0,i]=1}
data[,facratings:=NULL]
data <- cbind(data,dum)

# Add (lags of) rating upgrades and downgrades
data[,dratings0 := sign(get(idratings)-shift(get(idratings),n=1,type="lag")),by=get(idc)]
for (i in 1:lags_dratings){
  data[,paste0("dratings",i) := shift(dratings0,n=i,type="lag"),by=get(idc)]
}

if (cdums){
  # add country fixed effects
  coun <- dummy(idc,data)
  data <- cbind(data,coun)
}

if (asym) {
  # add extra columns for positive and negative changes
  for (choice in c("dyie","drat")) {
    cn <- colnames(data)
    i <- which(substr(cn,1,4)==choice)
    pos <- data[,..i]
    colnames(pos) <- paste(cn[i],"p",sep="_")
    neg <- data[,..i]
    colnames(neg) <- paste(cn[i],"n",sep="_")
    neg[neg>0] = 0
    pos[pos<0] = 0
    data <- cbind(data,pos,neg)
  }
}
data[,dratings0_fac := as.factor(dratings0)]

# Transform to data.frame for ease of calculations
data <- as.data.frame(data)
return(data)
}