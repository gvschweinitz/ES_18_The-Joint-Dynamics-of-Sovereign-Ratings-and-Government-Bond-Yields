IRF_mult <- function(data,tab_y,tab_r,std_y=0,fspecs=fspecs,empirical_changes=NULL,start_combis,downgrade = 2,conf_lvls=NULL,periods=120,diff_to_benchmark=TRUE){
  #DESCRIPTION:
  # Calculates IRFs based on bootstrap coefficient draws for a given downgrade scenario using IRF_yshock_bs.R
  # IRFs are calculated for different starting ratings
  #-------------------------------------------------------------------------------
  #USAGE
  # IRF_mult(data,tab_y,tab_r,std_y=0,fspecs=fspecs,empirical_changes=NULL,start_combis,downgrade = 2,conf_lvls=NULL,periods=120,diff_to_benchmark=TRUE)
  #-------------------------------------------------------------------------------
  #INPUT
  # data              data frame with data with the naming convention from data_prep
  # tab_y             bootstrap coefficients from yield equation
  # tab_r             bootstrap coefficients from rating equation
  # std_y             Standard deviation of yield changes. Needed to account for the influence of past changes due to asymmetric effects
  # fspecs            Structure controlling estimation setup
  # empirical_changes Result from extract_changes.R
  # start_combis      Matrix with columns "r" and "y given r". Gives (multiple) rating-yield combinations before the shock.
  # downgrade         Downgrade scenario.
  # conf_lvls         Confidence levels for IRFs. May also contain strings indicating standard deviations -> median +/- one standard deviation
  # periods           IRF horizon
  # diff_to_benchmark Boolean. FALSE: IRF relative to initial rating. TRUE: IRF relative to usual adjustment towards long-run equilibrium
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   start_combis      K different initial scenarios
  #   downgrade         downgrade scenario
  #   text              text description of downgrade scenario for each scenario
  #   IRF_mult_med_y    K x (periods+length(downgrade)+1) matrix of IRFs of yields
  #   IRF_mult_med_r    K x (periods+length(downgrade)+1) matrix of IRFs of ratings
  #   IRF_mult_sharedefault  K x (periods+length(downgrade)+1) matrix of default probability
  #   conf_lvls         Confidence levels for IRFs.
  #   IRF_mult_conf_y   K x (periods+length(downgrade)+1) x length(conf_lvl) array of IRFs of yields at all conf_lvls
  #   IRF_mult_conf_r   K x (periods+length(downgrade)+1) x length(conf_lvl) array of IRFs of ratings at all conf_lvls
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
  sdcheck <- c("+sd","+std","+SD","+STD","-sd","-std","-SD","-STD")
  
  minrat <- 6
  maxrat <- 24
  downgrade <- c(0,downgrade)
  start_p <- length(downgrade)
  pos <- which((start_combis[,"r"]>=(minrat+max(cumsum(downgrade)))) & (start_combis[,"r"]<=(maxrat+min(cumsum(downgrade)))))
#   start_combis <- start_combis[pos,]
  nscens <- length(pos)
  
  IRF_mult_med_y <- mat.or.vec(nscens,periods+start_p)
  IRF_mult_med_r <- mat.or.vec(nscens,periods+start_p)
  IRF_mult_sharedefault <- mat.or.vec(nscens,periods+start_p)+1
  if (!is.null(conf_lvls)){
    IRF_mult_conf_y <- array(0,dim=c(nscens,periods+start_p,length(conf_lvls)))
    IRF_mult_conf_r <- array(0,dim=c(nscens,periods+start_p,length(conf_lvls)))
  }
  
  for (n in 1:nscens){
    start <- NULL
    start$ratings <- start_combis[pos[n],"r"]-cumsum(downgrade)
    rat <- c(start$ratings[1],start$ratings[1:(start_p-1)])
    start$yields <- start_combis[match(rat,start_combis[,"r"]),"y given r"]
    print(paste("starting scenario",n))
    print(start)
    temp <- IRF_yshock_bs(data,tab_y,tab_r,std_y=std_y,fspecs=fspecs,empirical_changes=empirical_changes,start=start,periods=120,diff_to_benchmark=diff_to_benchmark)
    
    # consistent IRF closest to the median
    med_y <- apply(temp$IRF_y,2,quantile,0.5)
    utility <- rep(0,dim(temp$IRF_y)[1])
    for (i in 1:dim(temp$IRF_y)[1]){
      utility[i] <- sum((temp$IRF_y[i,]-med_y)^2)
    }
    pos_y <- which(utility == min(utility))[1]
    IRF_mult_med_y[n,] <- temp$IRF_y[pos_y,]
    
    med_r <- apply(temp$IRF_r,2,quantile,0.5)
    utility <- rep(0,dim(temp$IRF_r)[1])
    for (i in 1:dim(temp$IRF_r)[1]){
      utility[i] <- sum((temp$IRF_r[i,]-med_r)^2)
    }
    pos_r <- which(utility == min(utility))[1]
    IRF_mult_med_r[n,] <- temp$IRF_r[pos_r,]
    IRF_mult_sharedefault[n,] <- apply(temp$Default,2,mean)
    
    # confidence levels of IRF
    if (!is.null(conf_lvls)) {
      for (i in 1:length(conf_lvls)){
        if (!is.na(as.numeric(conf_lvls[i]))){
          IRF_mult_conf_y[n,,i] <- apply(temp$IRF_y,2,quantile,as.numeric(conf_lvls[i]))
          IRF_mult_conf_r[n,,i] <- apply(temp$IRF_r,2,quantile,as.numeric(conf_lvls[i]))
        }
        else{
          if (sum(grepl(conf_lvls[i],sdcheck))>0){
            if (!is.na(charmatch("+",conf_lvls[i]))) {
              IRF_mult_conf_y[n,,i] <- apply(temp$IRF_y,2,quantile,0.5) + apply(temp$IRF_y,2,sd)
              IRF_mult_conf_r[n,,i] <- apply(temp$IRF_r,2,quantile,0.5) + apply(temp$IRF_r,2,sd)
            }
            else {
              IRF_mult_conf_y[n,,i] <- apply(temp$IRF_y,2,quantile,0.5) - apply(temp$IRF_y,2,sd)
              IRF_mult_conf_r[n,,i] <- apply(temp$IRF_r,2,quantile,0.5) - apply(temp$IRF_r,2,sd)
            }
          }
        }
      }
    }
  }
  
  
  res <- NULL
  res$start_combis <- start_combis
  res$downgrade <- downgrade
  res$text <- paste(translate_rating(start_combis[pos,"r"]),translate_rating(start_combis[pos,"r"]-cumsum(downgrade[length(downgrade)])),sep=" -> ")
  res$IRF_mult_med_y <- IRF_mult_med_y
  res$IRF_mult_med_r <- IRF_mult_med_r
  res$IRF_mult_sharedefault <- IRF_mult_sharedefault
  if (!is.null(conf_lvls)){
    pos_sd <- grepl("sd",conf_lvls)
    if (sum(pos_sd)>0){conf_lvls <- conf_lvls[]}
    res$conf_lvls <- conf_lvls
    res$IRF_mult_conf_y <- IRF_mult_conf_y
    res$IRF_mult_conf_r <- IRF_mult_conf_r
  }
  return(res)
}