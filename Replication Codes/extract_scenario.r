extract_scenario <- function(month,year,country,periods = 6, data, res_yields, res_ratings,rat_var = "ratings_f_last", yield_var = "ryields") {
  #DESCRIPTION:
  # Extract real-world scenario (starting values) for simulation in IRF_scenario_par.R
  #-------------------------------------------------------------------------------
  #USAGE
  # extract_scenario(month,year,country,periods = 6, data, res_yields, res_ratings,rat_var = "ratings_f_last", yield_var = "ryields")
  #-------------------------------------------------------------------------------
  #INPUT
  # month             Month of the first observation of the real-world scenario
  # year              Year of the first observation of the real-world scenario
  # country           Country of the real-world scenario
  # periods           Length of the real-world scenario before simulation start
  # data              data frame with data with the naming convention from data_prep
  # res_yields        results from yield estimation
  # res_ratings       results from rating estimation
  # rat_var           name of rating variable
  # yield_var         name of yield variable
  #-------------------------------------------------------------------------------
  #OUTPUT
  # list with
  #   yield       yields of starting scenario
  #   ratings     ratings of starting scenario
  #   yield_X     data of starting scenario for yield estimation
  #   rating_X    data of starting scenario for rating estimation
  #   y_change    yield changes of starting scenario (LHS yield equation)
  #   r_change    rating changes of starting scenario (LHS rating equation)
  #   obs_yield   period_sim months of observed yields after initial scenario (for comparison with simulated IRF)
  #   ratings     period_sim months of observed ratings after initial scenario (for comparison with simulated IRF)
  #   period_sim  Number of periods for simulated IRF, including initial scenario
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
dateid = "dateid"
crossid = "country"
date_col <- as.Date(data[,dateid],"%m/%d/%Y")

if (is.character(country)) {country_col <- as.character(data[,crossid])}
if (is.numeric(country)) {country_col <- as.numeric(data[,crossid])}

T <- length(date_col)

month_vec <- as.numeric(format(date_col,"%m"))
year_vec <- as.numeric(format(date_col,"%Y"))

d <- which((month_vec==month) & (year_vec ==year))
c <- which(country_col == country)

i <- intersect(d,c)

i <- i:(i+periods-1)
print(i)
y_select <- res_yields$res_smooth$model[i,]
y_change <- y_select[,1]
y_select <- y_select[,-1]
r_select <- res_ratings$pre$X[i,]
r_change <- res_ratings$pre$Y[i]-2



result <- NULL
result$yield_X <- y_select
result$rating_X <- r_select
result$yield <- data[i,yield_var]
result$ratings <- data[i,rat_var]
result$r_change <- r_change
result$y_change <- y_change

pos_obs <- i[1]:max(c)
if (length(pos_obs)>120){pos_obs<-pos_obs[1:120]}
result$obs_yield <- data[pos_obs,yield_var]
result$obs_ratings <- data[pos_obs,rat_var]
result$period_sim <- length(pos_obs)
return(result)
}