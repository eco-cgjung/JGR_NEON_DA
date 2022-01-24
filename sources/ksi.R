ksi <- function(Nt,temp, moist, mscut, Q10) {
  ksi <- c()
  tmp <- c()
  #moisture <- c()
  for (i in 1:Nt) {
    sumtemp <- 0
    tmp <- temp[i]
    if (i > 10) {
      for (j in (i-9):i) {
        sumtemp <- sumtemp+temp[j]
      }
    tmp <- sumtemp/10
  }
  tmp <- Q10^((tmp-10)/10)
  moisture <- 1
  if (moist[i]<mscut) {
    moisture <- 1.0 - 5*(mscut-moist[i])
  }
  ksi[i] <- tmp*moisture
  }
  ksi
}