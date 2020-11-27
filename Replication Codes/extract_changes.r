extract_changes <- function(data, id= "crossid", rat = "ratings_f_last") {
  # DESCRIPTION:
  # collects empirical distribution of rating changes
  #-------------------------------------------------------------------------------
  #USAGE
  # extract_changes(data, id= "crossid", rat = "ratings_f_last")
  #-------------------------------------------------------------------------------
  #INPUT
  # data          data frame with data with the naming convention from data_prep
  # crossid       name of country column (N cross sections)
  # rat           name of rating information
  #-------------------------------------------------------------------------------
  #OUTPUT
  # matrix with
  #   (1) observed rating steps (up or down)
  #   (2) cumulative empirical distribution of upgrades
  #   (3) cumulative empirical distribution of downgrades
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  
ids <- unique(data[,id])

df <- NULL

for (n in ids) {
  d <- data[data[,id]==n,]
  df <- c(df,diff(d[,rat]))
}

df <- round(df*3)/3
df <- df[df!=0]

steps <- unique(abs(df))
steps <- sort(steps)

table <- mat.or.vec(length(steps),3)
table[,1] <- steps

for (n in 1:length(steps)) {
  table[n,2] <- sum(df==steps[n], na.rm = TRUE)
  table[n,3] <- sum(df==-steps[n], na.rm = TRUE)
}
table[,2] = cumsum(table[,2])/sum(table[,2])
table[,3] = cumsum(table[,3])/sum(table[,3])

return(table)
}