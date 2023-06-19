eweibull=function (probs, k, niter = 50) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (k <= 0) {
    stop("k must be greater than 0.")
  }
  e = rep(gamma(1 + 1/k), length(probs))
  i = 1
  while (i <= niter) {
    e = gamma(1 + 1/k) * ((2 * probs - 1) * (1-pgamma(e^k, 1 + 1/k)) + 1 - probs)/((2 * probs - 1) * exp(-e^k) + 
                                                                                     1 - probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  return(e)
}