einvgamma=function (probs, alpha, lambda, niter = 50) 
{
  if (min(probs) <= 0 || max(probs) > 1) {
    stop("only asymmetries between 0 (strictly) and 1 allowed.")
  }
  if (alpha <= 1 || lambda <= 0) {
    stop("alpha and beta must be strictly greater than 1 and 0, respectively.")
  }
  e = rep(lambda/(alpha-1), length(probs))
  i = 1
  while (i <= niter) {
    e = lambda/(alpha-1)*((2*probs-1)*pgamma(lambda/e,alpha-1)+1-probs)/((2*probs-1)*pgamma(lambda/e,alpha)+1-probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  return(e)
}
