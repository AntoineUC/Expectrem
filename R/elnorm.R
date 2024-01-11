elnorm=function (probs, mu = 0, sigma = 1, tol=1e-08,maxiter=100) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (sigma <= 0) {
    stop("sigma must be positive.")
  }
  e = rep(exp(mu+sigma^2/2), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<=maxiter) {
    e1 = exp(mu+sigma^2/2)*((2*probs-1)*(1-(pchisq(2 * ((log(e)-mu-sigma^2)/(sqrt(2)*sigma))^2, 1) * sign(((log(e)-mu-sigma^2)/(sqrt(2)*sigma)))))+2*(1-probs))/((2*probs-1)*(1-(pchisq(2 * ((log(e)-mu)/(sqrt(2)*sigma))^2, 1) * sign(((log(e)-mu)/(sqrt(2)*sigma)))))+2*(1-probs))
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
