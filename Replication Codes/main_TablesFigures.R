# This script creates all Tables and Figures in 
# "The joint dynamics of sovereign ratings and government bond yields"
# Makram El-Shagi and Gregor von Schweinitz, Journal of Banking and Finance (2018)

library("plyr")
wd <- getwd()
load(paste(wd,"/ratingfirst_asym.RData",sep=""))
source(paste(wd,"/estim_smooth_post_bs.R",sep=""))
do.col <- FALSE               # TRUE: colored figures
mai.default <- par("mai")     # Needed for resetting plot settings

#--------------------------------------------
#SUMMARY STATISTICS
IMFmatch <- matrix("",46,2)
IMFmatch[,1] <- as.character(unique(data$country))
IMFmatch[,2] <- "Advanced"
IMFmatch[c(1,5,7:9,18,19,24,25,29:31,39,43,44),2] <- "Developing"
IMFmatch[c(16,32,35,36),2] <- "Transition"
data$IMFclass <- IMFmatch[pmatch(data$country,IMFmatch[,1],duplicates.ok=TRUE),2]

# Table A2 and A3
sumstat_country <- ddply(data,.(country),summarize,
                         mindate=min(as.Date(dateid,format="%m/%d/%Y")),maxdate=max(as.Date(dateid,format="%m/%d/%Y")),
                         meany=mean(lyields,na.rm=TRUE),sdy=sd(lyields,na.rm=TRUE),miny=min(lyields,na.rm=TRUE),maxy=max(lyields,na.rm=TRUE),
                         meanr=mean(ratings_f_first,na.rm=TRUE),sdr=sd(ratings_f_first,na.rm=TRUE),minr=min(ratings_f_first,na.rm=TRUE),maxr=max(ratings_f_first,na.rm=TRUE),
                         maxdiff = max(diff(ratings_f_first,lag=12),na.rm=TRUE),mindiff = min(diff(ratings_f_first,lag=12),na.rm=TRUE),sddiff = sd(diff(ratings_f_first,lag=12),na.rm=TRUE))
sumstat_country$text <- paste(sumstat_country$mindate,"-",sumstat_country$maxdate)

# Table 2: Summary statistics by country group
sumstat_group <- data.frame(Indicator = "yields", IMFclass="total sample", mean=mean(data$lyields,na.rm=TRUE),sd=sd(data$lyields,na.rm=TRUE),min=min(data$lyields,na.rm=TRUE),max=max(data$lyields,na.rm=TRUE))
sumstat_group <- rbind(sumstat_group,
                       cbind(Indicator="yields",ddply(data,.(IMFclass),summarize,
                                                      mean=mean(lyields,na.rm=TRUE),sd=sd(lyields,na.rm=TRUE),min=min(lyields,na.rm=TRUE),max=max(lyields,na.rm=TRUE))))

sumstat_group <- rbind(sumstat_group,
                       cbind(Indicator = "ratings", IMFclass="total sample", mean=mean(data$ratings_f_first,na.rm=TRUE),sd=sd(data$ratings_f_first,na.rm=TRUE),min=min(data$ratings_f_first,na.rm=TRUE),max=max(data$ratings_f_first,na.rm=TRUE)))
sumstat_group <- rbind(sumstat_group,
                       cbind(Indicator="ratings",ddply(data,.(IMFclass),summarize,
                                                       mean=mean(ratings_f_first,na.rm=TRUE),sd=sd(ratings_f_first,na.rm=TRUE),min=min(ratings_f_first,na.rm=TRUE),max=max(ratings_f_first,na.rm=TRUE))))


##############
#SUMMARY PLOTS OF RATING AND YIELD DATA
# Figure 1: Density of yields
pdf("dens_yields.pdf",width=7,height=4)
d <- density(data$ryields,na.rm=TRUE)
plot(d$x,d$y,main="",xlab="Real yields (in per cent p.a.)",ylab="Density",type='l',lwd=2)
dev.off()

# Figure 2: Scatterplot of ratings and yields
pdf("ry_scat.pdf",width=7.5,height=4)
x <- translate_rating(c(6:24))
plot(-data$ratings_f_first,data$ryields,type="p",pch=20,main="",xlab="Ratings",ylab="Yields (in per cent p.a.)",xaxt='n')
abline(h=0)
axis(1,at=-c(6:24),x,cex.axis=1.1)
dev.off()

# Figure 3: Histogram of monthly ratings, by country group
pdf("hist_ratings_bygroup.pdf",width=7,height=4)
par(mai=rep(1.2, 4))
freqall <- ddply(data,.(IMFclass,lgratings),summarise,freq=length(lgratings))
freqall <- freqall[!is.na(freqall[,2]),]
rat <- seq(max(freqall[,2]),min(freqall[,2]),-1)
x <- translate_rating(rat)
freq <- matrix(0,length(rat),4)
freq[,1] <- rat
freq[match(freqall[freqall[,1]=="Advanced",2],freq),2] <- freqall[freqall[,1]=="Advanced",3]
freq[match(freqall[freqall[,1]=="Transition",2],freq),3] <- freqall[freqall[,1]=="Transition",3]
freq[match(freqall[freqall[,1]=="Developing",2],freq),4] <- freqall[freqall[,1]=="Developing",3]
colnames(freq) <- c("Level","Advanced","Transition","Developing")
barplot(height=t(freq[,2:4]),names.arg=x,legend.text=TRUE,args.legend=c(x="topright"),xlab="Ratings",ylab="Frequency",cex.axis=1.1,cex.lab=1.1)
par(mai=mai.default)
dev.off()

# Figure 4: Histogram of monthly rating changes, by country group
pdf("hist_rchanges_bygroup.pdf",width=7,height=4)
d <- ddply(data,.(country),transform,dratings = c(NA,diff(ratings_f_first,lag=1)))
d <- d[,c("IMFclass","dratings")]
d <- d[!is.na(d[,2]),]
d[,2] <- round(d[,2]*3)
d <- d[d[,2]!=0,]
freqall <- ddply(d,.(IMFclass,dratings),summarise,freq=length(dratings))
freqall <- freqall[!is.na(freqall[,2]),]
rat <- c(min(d[,2]):max(d[,2]))
freq <- matrix(0,length(rat),4)
freq[,1] <- rat
freq[match(freqall[freqall[,1]=="Advanced",2],freq),2] <- freqall[freqall[,1]=="Advanced",3]
freq[match(freqall[freqall[,1]=="Transition",2],freq),3] <- freqall[freqall[,1]=="Transition",3]
freq[match(freqall[freqall[,1]=="Developing",2],freq),4] <- freqall[freqall[,1]=="Developing",3]
freq[,1] <- freq[,1]/3
colnames(freq) <- c("Level","Advanced","Transition","Developing")
barplot(height=t(freq[,2:4]),names.arg=freq[,1],legend.text=TRUE,args.legend=c(x="topleft"),xlab="Rating Changes",ylab="Frequency",cex.axis=1.1)
dev.off()

# Figure 5: Event study, yield development before / after rating change
pdf("event_study.pdf",width=7,height=4)
Tdiff <- 12
d <- data[,c("country","dateid","yields_m","dratings0")]
countries <- unique(d[,"country"])
matpos <- matrix(NA,sum(d[,4]>0,na.rm=TRUE),2*Tdiff+1)
matneg <- matrix(NA,sum(d[,4]<0,na.rm=TRUE),2*Tdiff+1)
countpos <- 0
countneg <- 0
for (coun in countries){
  print(coun)
  temp <- d[d[,"country"]==coun,]
  T <- dim(temp)[1]
  pos <- which(temp[,4]!=0)
  for (i in pos){
    vec <- c(rep(NA,max(Tdiff-i+1,0)),temp[max(1,i-Tdiff):min(T,(i+Tdiff)),3],rep(NA,max(i+Tdiff-T,0)))
    vec <- vec/vec[Tdiff+1]*100
    if (temp[i,4]>0){
      countpos <- countpos+1
      matpos[countpos,] <- vec
    }
    else {
      countneg <- countneg+1
      matneg[countneg,] <- vec
    }
  }
}
matpos <- matpos[1:countpos,]
matneg <- matneg[1:countneg,]

plotpos <- apply(matpos,2,quantile,c(0.5),na.rm=TRUE)
plotneg <- apply(matneg,2,quantile,c(0.5),na.rm=TRUE)
plotvals <- rbind(plotpos,plotneg)
matplot(c(-Tdiff:Tdiff),t(plotvals),type='l',lty=c(2,3),lwd=c(2,2),col=1,xlab="Months before/after rating change",ylab="normalized yield",cex.axis=1.1)
legend("bottomleft",legend=c("yields around upgrade","yields around downgrade"),lty=c(2,3),lwd=c(2,2),col=1)
dev.off()

# Figure 6: idealized relations of ratings and yields
pdf("ideal_relation_rest.pdf",width=6,height=4)
par(mai=c(rep(0.8, 3),0.25))
x <- seq(-5,2,0.1)
y <- mat.or.vec(length(x),2)

y[,1] <- exp(x)
y[,2] <- 2*exp(x/2)-0.5
if (do.col){
  matplot(x,y,type='l',lwd=2,col=c(2,4),xaxt='n',yaxt='n',xlab="",ylab="")
}else{
  matplot(x,y,type='l',lwd=2,col=1,xaxt='n',yaxt='n',xlab="",ylab="")
}
title(xlab="Ratings",ylab="Yields",line=1,cex.lab=1)
xtext = c(-2.5,0.8)
ytext = c(0.5,3.2)
stext = c("good","bad")
text(xtext,ytext,stext)
par(mai=mai.default)
dev.off()

#--------------------------------------------
#LONG-RUN RELATIONSHIP
# Table 3 (second column of coeff_y, and coeff_r plus stats_bs)
coeff_y <- cbind(y_jointlambda$res_smooth$coefficients,apply(coeffs_bs1000$tab_y,2,quantile,0.5),colMeans(coeffs_bs1000$tab_y>0))
colnames(coeff_y) <-  c("mean_coeff","median_coeff","p>0")
coeff_r <- cbind(r_jointlambda$res_smooth$coefficients,apply(coeffs_bs1000$tab_r,2,quantile,0.5),colMeans(coeffs_bs1000$tab_r>0))
colnames(coeff_r) <-  c("mean_coeff","median_coeff","p>0")

y_jointlambda$post_bs <- estim_smooth_post_bs(y_jointlambda,c(coeff_y[,2]),oprob=FALSE)
r_jointlambda$post_bs <- estim_smooth_post_bs(r_jointlambda,c(coeff_r[,2]),oprob=TRUE)
stats_bs <- rbind(lambda=c(y_jointlambda$lambda,r_jointlambda$lambda),
                  coefficients = c(y_jointlambda$k,r_jointlambda$k),
                  LL_data=c(y_jointlambda$post_bs$LL_out$LL_data,r_jointlambda$post_bs$LL_out$LL_data),
                  LL_smooth=c(y_jointlambda$post_bs$LL_out$LL_smooth,r_jointlambda$post_bs$LL_out$LL_smooth),
                  R2=c(y_jointlambda$post_bs$R2,r_jointlambda$post_bs$R2),
                  R2adj=c(y_jointlambda$post_bs$R2adj,r_jointlambda$post_bs$R2adj),
                  BIC=c(y_jointlambda$post_bs$BIC,r_jointlambda$post_bs$BIC),
                  AIC=c(y_jointlambda$post_bs$AIC,r_jointlambda$post_bs$AIC))

# Table 4
time_taken_medy <- t(as.matrix(apply(time_taken_yshock,2,median)/12))

# Note: Figure 7 is plotted as part of main_estimation.R

#--------------------------------------------
#IRF PLOTS
# Figures 8 and 9: Median IRFs after two-notch downgrade from different starting values
pos <- c(1:5,11,14,17)
# 8 curves per plot
if (do.col){
  # differentiate by color
  reps <- ceiling(length(pos)/8)
  col <- rep(c(1:8),reps)
  type <- sort(rep(1:reps,8))
  pch <- rep(NA,length(pos))
  col.scen <- c(1,4,4,4,4,4)
  lty.scen <- c(1,1,2,2,3,3)
}else{
  print("b/w plots only work for 8 curves")
  col <- 1
  type <- rep(1:4,2)
  col.points <- 1:4
  pch <- c(1:4,rep(NA,4))
  point.pos <- seq(2,122,10)
  col.scen <- 1
  lty.scen <- c(6,1,2,2,3,3)
}

pdf("IRF_y_lim.pdf",width=10.5,height=7.5)
matplot(1:122,t(IRF_eqyields_down$IRF_mult_med_y[pos,]),xaxt="n",type="l",col=col,lty=type,lwd=2,xlab="Months",ylab="Yields",main="",cex.axis=1.1,cex.lab=1.2)
if (!do.col){for (k in col.points){points(point.pos,t(IRF_eqyields_down$IRF_mult_med_y[pos[k],point.pos]),pch=pch[k],lwd=2)}}
axis(1,at=seq(2,122,12),labels = seq(0,120,12),cex.axis=1.1)
legend("topright",legend=IRF_eqyields_down$text[pos],col=col,lty=type,pch = pch,lwd=2)
dev.off()

pdf("IRF_r_lim.pdf",width=10.5,height=7.5)
matplot(1:122,t(IRF_eqyields_down$IRF_mult_med_r[pos,]),xaxt="n",type="l",col=col,lty=type,lwd=2,xlab="Months",ylab="Ratings",main="",cex.axis=1.1,cex.lab=1.2)
if (!do.col){
  for (k in col.points){
    vals.all <- c(IRF_eqyields_down$IRF_mult_med_r[pos[k],])
    p.add <- which(diff(vals.all)!=0)
    vals.add <- apply(matrix(p.add),1,FUN=function(x){mean(vals.all[x:(x+1)])})
    points(c(point.pos,p.add+0.5),c(vals.all[point.pos],vals.add),pch=pch[k],lwd=2)
  }
}
axis(1,at=seq(2,122,12),labels = seq(0,120,12),cex.axis=1.1)
legend("topleft",legend=IRF_eqyields_down$text[pos],col=col,pch = pch,lty=type,lwd=2)
dev.off()

########
# Figures 10 and 11: Scenario plots
scen_names <- ls()[grep("scen_",ls())]
scen_names <- c("scen_ITA","scen_GRE_I")
for (i in 1:length(scen_names)){
  pdf(paste(scen_names[i],".pdf",sep=""),width=10.5,height=5.2)
  
  print_lvl <- c(0.05,0.95)
  scen <- get(scen_names[i])
  
  periods <- length(scen$scen$obs_yield)
  xtext <- seq(as.Date(paste(scen$year,"/",scen$month,"/1",sep="")),by="month",length.out=periods)
  pos_xtext <- which(format(xtext,"%m")=="01")
  text <- c("Observed development","IRF",
            paste(print_lvl[1],"conf without shocks"),paste(print_lvl[2],"conf without shocks"),
            paste(print_lvl[1],"conf with shocks"),paste(print_lvl[2],"conf with shocks"))
  text_exp <- paste("Scenario analysis for",scen$country,"in",xtext[scen$periods_start])
  pos_lvl <- charmatch(print_lvl,scen$conf_lvls)
  yields <- cbind(scen$scen$obs_yield,scen$IRF_med_y[1:periods],
                  scen$IRF_conf_med_y[1:periods,pos_lvl],
                  scen$IRF_conf_shock_y[1:periods,pos_lvl])
  ratings <- cbind(scen$scen$obs_ratings,scen$IRF_med_r[1:periods],
                   scen$IRF_conf_med_r[1:periods,pos_lvl],
                   scen$IRF_conf_shock_r[1:periods,pos_lvl])
  
  layout(matrix(c(1,2,3,3),ncol=2, byrow = TRUE), heights=c(4, 1.2))
  par(mai=c(1,1,1,0.3))
  matplot(x=c(1:periods),y=yields,type ="l",main="(a)",xaxt="n",xlab = "Months",ylab = "Yields",lwd=c(2,2,1,1,1,1),col = col.scen,lty = lty.scen,cex.axis=1.1,cex.lab=1.2)
  axis(1,at=pos_xtext,labels=format(xtext[pos_xtext],"%Y")) 
  abline(v=scen$periods_start)
  
  matplot(x=c(1:periods),y=ratings,type ="l",main="(b)",xaxt="n",xlab = "Months",ylab = "Ratings",lwd=c(2,2,1,1,1,1),col = col.scen,lty = lty.scen,cex.axis=1.1,cex.lab=1.2)
  axis(1,at=pos_xtext,labels=format(xtext[pos_xtext],"%Y"))
  abline(v=scen$periods_start)
  
  par(mai=c(0.5,0,0,0))
  plot.new()
  legend(x="center", ncol=2,legend=text,lwd=c(2,2,1,1,1,1),col = col.scen,lty = lty.scen)
  mtext(text_exp, side = 1, line = 0.5, outer = FALSE)
  par(mai=mai.default)
  dev.off()
}

##########
# Figure A1: IRF plots- median + conf for multiple rating shocks
posvec <- c(2,4,17)
print_lvl <- c(0.05,0.95)
pos_lvl <- charmatch(print_lvl,IRF_eqyields_down$conf_lvls)
# main <- c("(a)","(b)","(c)")

text <- c("Median IRF",paste(round(print_lvl[1]*100,0),"% confidence"),paste(round(print_lvl[2]*100,0),"% confidence"))
pdf("IRF_conf.pdf",width=10,height=10)
layout(matrix(c(c(1:6),7,7),ncol=2, byrow = TRUE), heights=c(3,3,3,1))
par(mai=c(0.7,0.8,0.5,0.3))
for (i in 1:length(posvec)){
  pos <- posvec[i]
  main <- paste(IRF_eqyields_down$text[pos],c(", Yields",", Ratings"))
  
  y_y <- cbind(IRF_eqyields_down$IRF_mult_med_y[pos,],IRF_eqyields_down$IRF_mult_conf_y[pos,,pos_lvl])
  matplot(x=c(1:122),y=y_y,xaxt="n",type ="l",main=main[1],xlab = "Months",ylab = "Yields",lwd=c(2,1,1),col=1,lty=c(1,2,2),cex.lab=1.2)
  axis(1,at=seq(2,122,12),labels = seq(0,120,12))
  abline(a=0,b=0)
  
  y_r <- cbind(IRF_eqyields_down$IRF_mult_med_r[pos,],IRF_eqyields_down$IRF_mult_conf_r[pos,,pos_lvl])
  matplot(x=c(1:122),y=y_r,xaxt="n",type ="l",main=main[2],xlab = "Months",ylab = "Ratings",lwd=c(2,1,1),col=1,lty=c(1,2,2),ylim=c(-3,1),cex.lab=1.2)
  axis(1,at=seq(2,122,12),labels = seq(0,120,12))
  abline(a=0,b=0)
}
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=3,legend=text,lwd=c(2,1,1),col=1,lty=c(1,2,2),cex=1.2)
par(mai=mai.default)
dev.off()

##########
# Figures A2 and A3: IRF plots different shock scenarios
pos <- c(1:5,11,14,17)
for (k_scen in 1:2){
  if (k_scen==1){
    pdf("IRF_staggered.pdf",width=11,height=6.2)
    temp <- IRF_eqyields_staggered
  }else{
    pdf("IRF_up.pdf",width=11,height=6.2)
    temp <- IRF_eqyields_up
  }
  layout(matrix(c(1,2,3,3),ncol=2, byrow = TRUE), heights=c(5, 1.2))
  par(mai = c(0.8,0.8,1,0.2))
  T <- dim(temp$IRF_mult_med_y)[2]
  T.point <- seq(T%%12,T,12) 
  matplot(x=c(1:T),y=t(temp$IRF_mult_med_y[pos,]),xaxt="n",type ="l",main="IRF yields",xlab = "Months",ylab = "Yields",lwd=1,col=col,lty=type,cex.lab=1.2)
  if (!do.col){for (k in col.points){points(T.point,t(temp$IRF_mult_med_y[pos[k],T.point]),pch=pch[k])}}
  axis(1,at=seq(T%%12,T,12),labels = seq(0,T-T%%12,12),cex.axis=0.8)
  abline(a=0,b=0)
  
  matplot(x=c(1:T),y=t(temp$IRF_mult_med_r[pos,]),xaxt="n",type ="l",main="IRF ratings",xlab = "Months",ylab = "Ratings",lwd=1,col=col,lty=type,cex.lab=1.2)
  if (!do.col){
    for (k in col.points){
      vals.all <- c(temp$IRF_mult_med_r[pos[k],])
      p.add <- which(diff(vals.all)!=0)
      vals.add <- apply(matrix(p.add),1,FUN=function(x){mean(vals.all[x:(x+1)])})
      points(c(T.point,p.add+0.5),c(vals.all[T.point],vals.add),pch=pch[k],lwd=1)
    }
  }
  axis(1,at=seq(T%%12,T,12),labels = seq(0,T-T%%12,12),cex.axis=0.8)
  abline(a=0,b=0)
  
  par(mai=c(0.5,0,0,0))
  plot.new()
  legend(x="center",ncol=4,legend=paste(temp$text[pos],"   "),col=col,pch = pch,lty=type,lwd=1,cex=0.8)
  par(mai=mai.default)
  dev.off()
  rm(temp,T,T.point)
}

###################
# Figures A4 and A5: IRFs to robustness check with median ratings
rm(list=ls())
wd <- getwd()
load(paste(wd,"/medratingfirst_asym.RData",sep=""))

do.col <- FALSE
# 8 curves per plot
pos <- c(1:5,11,14,17)
if (do.col){
  # differentiate by color
  col <- c(1:8)
  type <- 1
  pch <- rep(NA,8)
}else{
  col <- 1
  type <- rep(1:4,2)
  col.points <- 1:4
  pch <- c(1:4,rep(NA,4))
  point.pos <- seq(2,122,10)
}

pdf("IRF_y_lim_med.pdf",width=10.5,height=7.5)
matplot(1:122,t(IRF_eqyields_down$IRF_mult_med_y[pos,]),xaxt="n",type="l",col=col,lty=type,lwd=2,xlab="Months",ylab="Yields",main="",cex.axis=1.1,cex.lab=1.2)
if (!do.col){for (k in col.points){points(point.pos,t(IRF_eqyields_down$IRF_mult_med_y[pos[k],point.pos]),pch=pch[k],lwd=2)}}
axis(1,at=seq(2,122,12),labels = seq(0,120,12),cex.axis=1.1)
legend("topright",legend=IRF_eqyields_down$text[pos],col=col,lty=type,pch = pch,lwd=2)
dev.off()

pdf("IRF_r_lim_med.pdf",width=10.5,height=7.5)
matplot(1:122,t(IRF_eqyields_down$IRF_mult_med_r[pos,]),xaxt="n",type="l",col=col,lty=type,lwd=2,xlab="Months",ylab="Ratings",main="",cex.axis=1.1,cex.lab=1.2)
if (!do.col){
  for (k in col.points){
    vals.all <- c(IRF_eqyields_down$IRF_mult_med_r[pos[k],])
    p.add <- which(diff(vals.all)!=0)
    vals.add <- apply(matrix(p.add),1,FUN=function(x){mean(vals.all[x:(x+1)])})
    points(c(point.pos,p.add+0.5),c(vals.all[point.pos],vals.add),pch=pch[k],lwd=2)
  }
}
axis(1,at=seq(2,122,12),labels = seq(0,120,12),cex.axis=1.1)
legend("topleft",legend=IRF_eqyields_down$text[pos],col=col,pch = pch,lty=type,lwd=2)
dev.off()