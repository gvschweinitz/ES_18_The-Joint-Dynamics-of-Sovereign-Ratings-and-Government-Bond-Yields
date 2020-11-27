translate_rating <- function(rat) {
  #DESCRIPTION
  # Translates numeric rating classes to descriptions.
  #-------------------------------------------------------------------------------
  #AUTHORS: Makram El-Shagi and Gregor von Schweinitz
  rat[rat<6] = 6
  rat[rat>24] = 24
  RAT <- c("AAA","AA+","AA","AA-","A+","A","A-","BBB+","BBB","BBB-","BB+","BB","BB-","B+","B","B-","CCC+","CCC","<CCC")
  return(RAT[25-rat])
}
  