einvgamma=function (probs, alpha, lambda, tol=1e-08, maxiter=100) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 (strictly) and 1 allowed.")
  }
  if (alpha <= 1 || lambda <= 0) {
    stop("alpha and beta must be strictly greater than 1 and 0, respectively.")
  }
  e = rep(lambda/(alpha-1), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<=maxiter) {
    e1 = lambda/(alpha-1)*((2*probs[which(probs*(1-probs) != 0)]-1)*pgamma(lambda/e[which(probs*(1-probs) != 0)],alpha-1)+1-probs[which(probs*(1-probs) != 0)])/((2*probs[which(probs*(1-probs) != 0)]-1)*pgamma(lambda/e[which(probs*(1-probs) != 0)],alpha)+1-probs[which(probs*(1-probs) != 0)])
    gap=max(abs(e1-e[which(probs*(1-probs) != 0)]),na.rm=T)
    e[which(probs*(1-probs) != 0)]=e1
    i=i+1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
