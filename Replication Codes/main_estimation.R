# This script runs all regressions for
# "The joint dynamics of sovereign ratings and government bond yields"
# Makram El-Shagi and Gregor von Schweinitz, Journal of Banking and Finance (2018)

#loading functions:
# IMPORTANT!! Adapt data_prep such that the right yields/ratings/... are used
wd <- getwd()
source(paste(wd, '/data_prep.R',sep=""))
source(paste(wd, '/run_BIC.R',sep=""))
source(paste(wd, '/cuttolags.R',sep=""))
source(paste(wd, '/eq_generator_bylag.R',sep=""))
source(paste(wd, '/eq_generator.R',sep=""))
source(paste(wd, '/IC_all_dyn.R',sep=""))
source(paste(wd, '/IC_all_dyn_joint.R',sep=""))
source(paste(wd, '/IC_all_dyn_coint.R',sep=""))
source(paste(wd, '/estim_smooth.R',sep=""))
source(paste(wd, '/estim_smooth_pre.R',sep=""))
source(paste(wd, '/estim_smooth_main.R',sep=""))
source(paste(wd, '/estim_smooth_post.R',sep=""))
source(paste(wd, '/LL_oprob.r',sep=""))
source(paste(wd, '/gr_oprob.r',sep=""))
source(paste(wd, '/hess_oprob.r',sep=""))
source(paste(wd, '/estim_coint.R',sep=""))
source(paste(wd, '/estim_coint_pre.R',sep=""))
source(paste(wd, '/estim_coint_main.R',sep=""))
source(paste(wd, '/estim_coint_post.R',sep=""))
source(paste(wd, '/LL_coint.R',sep=""))
source(paste(wd, '/extract_changes.r',sep=""))
source(paste(wd, '/bootstrap_sample.r',sep=""))
source(paste(wd, '/bootstrap_par.r',sep=""))
source(paste(wd, '/build_scen.R',sep=""))
source(paste(wd, '/extract_scenBS.r',sep=""))
source(paste(wd, '/simulation.r',sep=""))
source(paste(wd, '/analyze_sample.r',sep=""))
source(paste(wd, '/plot_longrun_bs.R',sep=""))
source(paste(wd, '/translate_rating.R',sep=""))
source(paste(wd, '/compare_models.R',sep=""))
source(paste(wd, '/IRF_mult.R',sep=""))
source(paste(wd, '/IRF_yshock_bs.R',sep=""))
source(paste(wd, '/sim_recovery_time_par.R',sep=""))
source(paste(wd, '/extract_scenario.r',sep=""))
source(paste(wd, '/IRF_scenario_par.R',sep=""))

# Option setting
fspecs <- NULL
fspecs$idc = "country"
fspecs$idobs = "dateid"
fspecs$idratings <- "ratings_f_first"          # Robustness: "ratings_med_first"
fspecs$idyields <- "yields_m"
fspecs$realyields <- TRUE
fspecs$iddeflator <- "cpi"
fspecs$idyields_final <- "ryields"
fspecs$yield_first <- FALSE
fspecs$na.rm <- TRUE
fspecs$cdums <- FALSE                          # Robustness: TRUE
fspecs$asym <- TRUE
fspecs$p_yields <- 1                   
fspecs$real_exp <- FALSE
fspecs$contemp <- TRUE

general <- "ratingfirst_asym"          #Baseline
# general <- "medratingfirst_asym"     #Robustness

# ------------------------------------------------------------------------------------
# DATA PREPARATION
data <- data_prep("data_monthly_Jan2016",12,12,fspecs = fspecs)

filenme <- paste(general,".RData",sep="")
longrun_plot <- paste("longrun_",general,".pdf",sep="")


# ------------------------------------------------------------------------------------
# LAG SELECTION by BIC in separate estimations (section 3.1.4)
# Lags in yield equation, with rating smoothing
optimcol <- 3
fspecs$factors <- c(7:max(data$gratings))
y_BIC <- run_BIC(data,maxlags=6,fspecs=fspecs,optimcol=optimcol)
fspecs$optlag_y <- which(y_BIC[,optimcol]==min(y_BIC[,optimcol]))
fspecs$ylags_y <- 1:fspecs$optlag_y
if (fspecs$yield_first){fspecs$ylags_r <- 1:(fspecs$optlag_y+1)}
if (!fspecs$yield_first){fspecs$ylags_r <- 0:fspecs$optlag_y}
eq_yields <- eq_generator_bylag(data,fspecs=fspecs,oprob=FALSE)

# Lags in rating equation, with full rating smoothing
r_BIC <- run_BIC(data,maxlags=6,fspecs=fspecs,oprob=TRUE,optimlambda=FALSE,optimcol=optimcol)
fspecs$optlag_r <- which(r_BIC[,optimcol]==min(r_BIC[,optimcol]))
fspecs$optlag_r <- 2  # ADJUSTMENT BECAUSE BIC is pretty equal for two and four lags
if (fspecs$yield_first){fspecs$rlags_y <- 0:fspecs$optlag_r}
if (!fspecs$yield_first){fspecs$rlags_y <- 1:(fspecs$optlag_r+1)}
fspecs$rlags_r <- 1:fspecs$optlag_r
eq_ratings <- eq_generator_bylag(data,fspecs=fspecs,oprob=TRUE)

data_all <- data
data <- cuttolags(data,max(fspecs$ylags_y,fspecs$rlags_y),max(fspecs$ylags_r,fspecs$rlags_r))
save.image(filenme)

# ------------------------------------------------------------------------------------
# LAMBDA SELECTION (optimal smoothing parameters) for joint estimations
# Variant 1: joint estimation, but separate long-run relationships
jointlambda_dyn <- IC_all_dyn_joint(eq_yields,eq_ratings,data,max_diff=0.01)
ljoint <- jointlambda_dyn[which(jointlambda_dyn[,4]==min(jointlambda_dyn[,4])),1]
save.image(filenme)
# Variant 2: joint estimation, but common long-run relationships (cointegration)


lcoint <- coint_BIC[which(coint_BIC[,3]==min(coint_BIC[,3])),1]
save.image(filenme)

# ------------------------------------------------------------------------------------
# MODEL COMPARISON: estimate models (separate versus common long-run relationship) with different candidates for smoothing parameter
ry_coint <- estim_coint(lcoint,data,fspecs$optlag_y,fspecs$optlag_r,cdums=fspecs$cdums,rasym=fspecs$asym,yasym=fspecs$asym)
ry_coint_alt <- estim_coint(ljoint,data,fspecs$optlag_y,fspecs$optlag_r,cdums=fspecs$cdums,rasym=fspecs$asym,yasym=fspecs$asym)

y_jointlambda <- estim_smooth(eq_yields,data,ljoint)
y_baseline <- estim_smooth(eq_yields,data,y_BIC[fspecs$optlag_y,1])
y_jointlambda_alt <- estim_smooth(eq_yields,data,lcoint)

r_jointlambda <- estim_smooth(eq_ratings,data,ljoint,oprob=TRUE)
r_baseline <- estim_smooth(eq_ratings,data,r_BIC[fspecs$optlag_r,1],oprob=TRUE)
r_jointlambda_alt <- estim_smooth(eq_ratings,data,lcoint,oprob=TRUE)

compare_models(y_jointlambda,r_jointlambda,ry_coint)
compare_models(y_jointlambda,r_jointlambda,ry_coint_alt)
compare_models(y_jointlambda_alt,r_jointlambda_alt,ry_coint)
save.image(filenme)

# ------------------------------------------------------------------------------------
# BOOTSTRAP ESTIMATION
# (1) Remove outliers
# Criterion: positive errors in yield equation with a probability of below 10^-8 under normally distributed errors
std_eps_full <- sqrt(y_jointlambda$post$var_data)
pos_outlier <- which(abs(y_jointlambda$res_smooth$residuals)>-qnorm(10^-8,sd=std_eps_full))
outlier <- cbind(data[pos_outlier,c("dateid","country","ryields","ratings_f_first")],eps_y=y_jointlambda$res_smooth$residuals[pos_outlier])
data_full <- data
data <- data_full[-pos_outlier,]
std_y=sd(data[,"dyields0"],na.rm=TRUE)
save.image(filenme)

# (2) Run bootstrap
sample_info_jointlambda <- analyze_sample(data,y_jointlambda)
empirical_changes <- extract_changes(data)
coeffs_bs1000 <- bootstrap_par(data,y_jointlambda,r_jointlambda,eq_yields,eq_ratings,fspecs=fspecs, sample_info=sample_info_jointlambda)
# Increase number of draws to be sure
# coeffs_bs5000 <- bootstrap_par(data,y_jointlambda,r_jointlambda,eq_yields,eq_ratings,fspecs=fspecs, sample_info=sample_info_jointlambda, iterations=5000)
# test_bs <- list(test_y = rep(NA,dim(coeffs_bs1000$tab_y)[2]),test_r = rep(NA,dim(coeffs_bs1000$tab_r)[2]))
# for (k in 1:dim(coeffs_bs1000$tab_y)[2]){
#   test_bs$test_y[k] <- ks.test(coeffs_bs1000$tab_y[,k],coeffs_bs5000$tab_y[,k])$p
# }
# for (k in 1:dim(coeffs_bs1000$tab_r)[2]){
#   test_bs$test_r[k] <- ks.test(coeffs_bs1000$tab_r[,k],coeffs_bs5000$tab_r[,k])$p
# }
# rm(k)
save.image(filenme)


# -----------------------------------------------------------------------------------
#RESULTS
# Plot Figure 7 and calculate long-run relationship
eq_plot <- plot_longrun_bs(coeffs_bs1000$tab_y,coeffs_bs1000$tab_r,coefficients_ML_y=NULL,coefficients_ML_r=NULL,std_y=std_y,savenme=longrun_plot)
save.image(filenme)

# Calculate IRFs for different downgrade scenarios (basis for IRF plots)
std_eps <- sd(y_jointlambda$res_smooth$residuals[1:length(y_jointlambda$pre$Y)])
# Case 1: rating downgrade by two notches
IRF_eqyields_down <- IRF_mult(data,coeffs_bs1000$tab_y,coeffs_bs1000$tab_r,std_y = std_eps,fspecs=fspecs,empirical_changes,eq_plot$yfair,conf_lvls=c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99,"-sd","+sd"),diff_to_benchmark=TRUE)
save.image(filenme)
# Case 2: rating downgrade by two notches, but staggered in two steps
IRF_eqyields_staggered <- IRF_mult(data,coeffs_bs1000$tab_y,coeffs_bs1000$tab_r,std_y = std_eps,fspecs=fspecs,empirical_changes,eq_plot$yfair,downgrade=c(1,1),conf_lvls=c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99,"-sd","+sd"),diff_to_benchmark=TRUE)
save.image(filenme)
# Case 3: rating upgrade by two notches
IRF_eqyields_up <- IRF_mult(data,coeffs_bs1000$tab_y,coeffs_bs1000$tab_r,std_y = std_eps,fspecs=fspecs,empirical_changes,eq_plot$yfair,downgrade=-2,conf_lvls=c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99,"-sd","+sd"),diff_to_benchmark=TRUE)
save.image(filenme)

# Estimate recovery times (basis for Table 4)
# time_taken_determ <- sim_recovery_time_par(y_jointlambda,r_jointlambda,fspecs=fspecs, start=c(rating=6,yield=eq_plot$yfair[1,2]),empirical_changes=empirical_changes, periods=6000,determ=TRUE)
time_taken_yshock <- sim_recovery_time_par(y_jointlambda,r_jointlambda,fspecs=fspecs, start=c(rating=6,yield=eq_plot$yfair[1,2]),empirical_changes=empirical_changes, periods=6000,determ=FALSE)
save.image(filenme)

# Estimate IRFs for real-world scenarios
scen_ITA <- IRF_scenario_par(5,1992,country="Italy",periods_start=5,data=data,fspecs=fspecs,res_yields=y_jointlambda,res_ratings=r_jointlambda,coeffs_bs=coeffs_bs1000,empirical_changes=empirical_changes,conf_lvls=c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99))
save.image(filenme)
scen_GRE_I <- IRF_scenario_par(9,2009,country="Greece",periods_start=11,data=data,fspecs=fspecs,res_yields=y_jointlambda,res_ratings=r_jointlambda,coeffs_bs=coeffs_bs1000,empirical_changes=empirical_changes,conf_lvls=c(0.01,0.05,0.1,0.2,0.5,0.8,0.9,0.95,0.99))
save.image(filenme)