elnorm=function (probs, mu = 0, sigma = 1, niter = 50) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (sigma <= 0) {
    stop("sigma must be positive.")
  }
  e = rep(exp(mu+sigma^2/2), length(probs))
  i = 1
  while (i <= niter) {
    e = exp(mu+sigma^2/2)*((2*probs-1)*(1-(pchisq(2 * ((log(e)-mu-sigma^2)/(sqrt(2)*sigma))^2, 1) * sign(((log(e)-mu-sigma^2)/(sqrt(2)*sigma)))))+2*(1-probs))/((2*probs-1)*(1-(pchisq(2 * ((log(e)-mu)/(sqrt(2)*sigma))^2, 1) * sign(((log(e)-mu)/(sqrt(2)*sigma)))))+2*(1-probs))
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  return(e)
}