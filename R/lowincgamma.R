lowincgamma=function (x, a) 
{
  if (!is.numeric(a) || !is.numeric(x)) 
    stop("All arguments must be real numbers.")
  
  xam <- -x + a * log(x)
  if (max(abs(xam)) > 700 || abs(a) > 170) {
    warning("Arguments 'x' and/or 'a' are too large.")
    return(NA)
  }
  gin=rep(0,length(x))
  
  s <- 1/a
  r <- s
  for (k in 1:60) {
    r <- r * x[which(x<=1+a)]/(a + k)
    s <- s + r
    if (max(abs(r/s)) < 1e-15) 
      break
  }
  gin[which(x<=1+a)] <- exp(xam[which(x<=1+a)]) * s
  
  t0 <- 0
  for (k in 60:1) {
    t0 <- (k - a)/(1 + k/(x[which(x>1+a)] + t0))
  }
  gim <- exp(xam[which(x>1+a)])/(x[which(x>1+a)] + t0)
  ga <- gamma(a)
  gin[which(x>1+a)] <- ga - gim
  
  return(Re(gin))
}