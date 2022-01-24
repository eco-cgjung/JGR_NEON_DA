##functions
roll.mean.temp <- function(temp) {
  roll.mean.temp <- rollmean(temp, 10)
}

GDD.func <- function(temp) {
  GDD <- cumsum(pmax((temp - 5), 0))
  # GDD <- c(rep(0,9),cumsum(pmax((temp-3), 0)))
}

FDD.func <- function(temp) {
  # FDD <- c(rep(22,9),cumsum(pmin((temp-22), 0)))
  FDD <- cumsum(pmin((temp-22), 0))
}

##switch
fun.sw2 <- function(GDD, roll.mean.T, T.fall, date.limit, offset) { #date.limit <- 200
  DOY.full <- length(GDD)
  sw <- ifelse((GDD <= GDD0), 1, offset)
  sw.2 <- ifelse(roll.mean.T < T.fall, 1, offset)
  sw[date.limit:DOY.full] <- sw.2[date.limit:DOY.full]
  sw <- c(rep(1,9), sw)
  sw
}