plot_longrun_bs <- function(coefficients_y,coefficients_r,coefficients_ML_y=NULL, coefficients_ML_r=NULL,std_y=0, savenme = NULL, conf_lev=NULL){
  #DESCRIPTION:
  # Plots long-run relationships of yield and rating equation, together with long-run equilibrium
  # Figure 7 in the paper
  #-------------------------------------------------------------------------------
  #USAGE
  # plot_longrun_bs(coefficients_y,coefficients_r,coefficients_ML_y=NULL, coefficients_ML_r=NULL,std_y=0, savenme = NULL, conf_lev=NULL)
  #-------------------------------------------------------------------------------
  #INPUT
  # coefficients_y      bootstrap coefficients from yield equation
  # coefficients_r      bootstrap coefficients from rating equation
  # coefficients_ML_y   ML estimates of coefficients from yield equation
  # coefficients_ML_r   ML estimates of coefficients from rating equation
  # std_y               Standard deviation of yield changes. Needed to account for the influence of past changes due to asymmetric effects
  # savenme             plot name
  # conf_lvl            Confidence bands to be plotted
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   yfair         matrix of fair yields for every rating in the yield and rating equation (median plots)
  #   pvals         matrix of p-values at every rating for positive yields in both equations, rating equation above yield equation and pressure away from equilibrium
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  save = ifelse(is.null(savenme),FALSE,TRUE)
  if(is.null(conf_lev)){conf_lev <- c(0.1,0.9)}
  yc_names <- colnames(coefficients_y)
  rc_names <- colnames(coefficients_r)
  N_ratings <- sum(grepl("facratings",rc_names))
  facs <- NULL
  fac_str <- strsplit(rc_names[grepl("facratings",rc_names)],"facratings")
  for (k in 1:N_ratings){
    facs <- c(facs,as.numeric(fac_str[[k]][2]))
  }
  
  posy_y <- grepl("lyields",yc_names)
  posy_r <- grepl("lyields",rc_names)
  posc_y <- grepl("const",yc_names)
  posup_r <- (length(rc_names)-1)
  poslo_r <- length(rc_names)
  posf_y <- grepl("facratings",yc_names)
  posf_r <- grepl("facratings",rc_names)
  poslag_y <- grep("dyields",yc_names)
  poslag_r <- grep("dyields",rc_names)
  poslag_yp <- poslag_y[grep("_p",yc_names[poslag_y])]
  poslag_yn <- poslag_y[grep("_n",yc_names[poslag_y])]
  poslag_rp <- poslag_r[grep("_p",rc_names[poslag_r])]
  poslag_rn <- poslag_r[grep("_n",rc_names[poslag_r])]
  
  Ny <- sum(posf_y)+1
  Nr <- sum(posf_r)+1
  Dy <- mat.or.vec(Ny,length(yc_names))
  Dy[,posc_y] <- 1
  Dy[2:Ny,posf_y] <- 1
  Dy[2:Ny,posf_y] <- Dy[2:Ny,posf_y] - upper.tri(Dy[2:Ny,posf_y])
  Dy[,poslag_yp] <- 0.5*std_y
  Dy[,poslag_yn] <- -0.5*std_y
  
  Dr <- mat.or.vec(sum(posf_r)+1,length(rc_names))
  Dr[,posup_r] <- -0.5
  Dr[,poslo_r] <- -0.5
  Dr[2:Nr,posf_r] <- 1
  Dr[2:Nr,posf_r] <- Dr[2:Nr,posf_r] - upper.tri(Dr[2:Nr,posf_r])
  Dr[,poslag_rp] <- 0.5*std_y
  Dr[,poslag_rn] <- -0.5*std_y
  
  colnames(Dy) <- yc_names
  colnames(Dr) <- rc_names
  
  yvals <- mat.or.vec(dim(coefficients_y)[1],Ny)
  rvals <- mat.or.vec(dim(coefficients_r)[1],Nr)
  for(t in 1:dim(coefficients_y)[1]){
    yvals[t,] <- (-1/coefficients_y[t,posy_y]) * Dy %*% c(coefficients_y[t,])
    rvals[t,] <- (-1/coefficients_r[t,posy_r]) * Dr %*% c(coefficients_r[t,])
  }
  yplot <- apply(yvals,2,quantile,conf_lev)
  rplot <- apply(rvals,2,quantile,conf_lev)
  if (is.null(coefficients_ML_y)){ymean <- apply(yvals,2,quantile,0.5)}
  else {ymean <- c((-1/coefficients_ML_y[posy_y]) * Dy %*% c(coefficients_ML_y))}
  if (is.null(coefficients_ML_r)){rmean <- apply(rvals,2,quantile,0.5)}
  else {rmean <- c((-1/coefficients_ML_r[posy_r]) * Dr %*% c(coefficients_ML_r))}
  if (save) {
    splits <- strsplit(savenme,".",fixed=TRUE)[[1]]
    if (splits[length(splits)]=="tiff"){
      tiff(savenme,width=8.5,height=6,units="in",compression="none",res=300)
    }else if (splits[length(splits)]=="eps"){
      setEPS()
      postscript(savenme,width=8.5,height=6)
    }else{
      pdf(savenme,width=8.5,height=6)
    }
  }
  
  y <- t(rbind(ymean,rmean,yplot,rplot))
  
  y_int <- approx(c(6,facs),ymean,n=1000)
  r_int <- approx(c(6,facs),rmean,n=1000)
  diff <- abs(y_int$y-r_int$y)
  plotpos <- c(y_int$x[diff==min(diff)],y_int$y[diff==min(diff)])
  text_lab <- paste("Rating=",translate_rating(round(plotpos[1])),"; Yield=",round(plotpos[2],2),"%",sep="")
  x <- translate_rating(c(6,facs))
  matplot(x=-c(min(facs)-1,facs),y=y,type ="l",main="",xlab = "Ratings",ylab = "Yields (%)",lwd=c(2,2,1,1,1,1),xaxt='n',col=1,lty=c(1,2,1,1,2,2))
  axis(1,at=-c(min(facs)-1,facs),labels=x,cex.axis=0.7)
  legend("topleft",legend=c("long run relation (yield equation)","long run relation (rating equation)"),lty=c(1,2),col=c(1,1),lwd=2)
  text(-plotpos[1], plotpos[2]+20, labels=text_lab, cex= 0.9, pos=3) 

  yfair <- rbind(c(6,facs),apply(yvals,2,quantile,0.5),apply(rvals,2,quantile,0.5))
  rownames(yfair) <- c("r","y given r","r given y")
  colnames(yfair) <- x
  pvals <- 1-rbind(colMeans(yvals>0),colMeans(rvals>0),colMeans(rvals>yvals),
                   colMeans(((rvals>yvals)&(rvals>0))|((rvals<yvals)&(rvals<0))))
  rownames(pvals) <- c("p yields","p ratings","p Difference","p pressure")
  colnames(pvals) <- x
  print(pvals)
  if (save) {dev.off()}
  
  out <- NULL
  out$yfair <- t(yfair)
  out$pvals <- t(pvals)
  return(out)
}