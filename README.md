# ES_18_The Joint Dynamics of Sovereign Ratings and Government Bond Yields
 Replication Code and Data for El-Shagi and von Schweinitz (JBF, 2018)


#################################
The files in these folders can be used to replicate the paper
"The joint dynamics of sovereign ratings and government bond yields"
by Makram El-Shagi and Gregor von Schweinitz

#################################
Please cite as 
El-Shagi, M., & von Schweinitz, G. (2018). The joint dynamics of sovereign ratings and government bond yields. Journal of Banking & Finance, 97, 198-218.

In case of questions and errors, please write an email to
gsz@iwh-halle.de

#################################
The folder "Rating Data" contains the raw rating information obtained from www.countryeconomy.com. 
This information can be used to construct further rating statistics

The folder "Replication codes" contains all codes to run the codes in the paper. To replicate, run
- main_estimation.R: runs all estimations
- main_TablesFigures.R: plots all Figures, collects all information for Tables.

Among others, the folder contains codes
- run_BIC.R/coint_BIC.R: estimate the optimal degree of smoothing
- estimate_smooth.R: estimate the linear yield regression and the ordered probit rating equation separately
- estimate_coint.R: estimate a cointegration model of both equations
- bootstrap_par.R: wild bootstrap to estimate uncertainty around estimated coefficients
- IRF_yshock_bs.R: simulate IRFs accounting for asymmetric coefficients and long-run convergence