generate.pars <- function(p_op,pmin,pmax,d) {
  while(TRUE){
    rand <- runif(length(pmin))
    p_new <- p_op+(rand-0.5)*(pmax-pmin)/d
    if  (Reduce("&", p_new>pmin&p_new<pmax))
      break
  }
  p_new
}

generate.pars.cov <- function(p_op, covars) {
  while(TRUE){
    p_new <- rmvn(1, mu = p_op, sigma = covars,ncores = 12)
    if  (Reduce("&", p_new>pmin&p_new<pmax))
      break
  }
  p_new
}