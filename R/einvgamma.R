einvgamma=function (probs, alpha, lambda, tol=1e-08, maxiter=100) 
{
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 (strictly) and 1 allowed.")
  }
  if (alpha <= 1 || lambda <= 0) {
    stop("alpha and beta must be strictly greater than 1 and 0, respectively.")
  }
  e = rep(lambda/(alpha-1), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<=maxiter) {
    e1 = lambda/(alpha-1)*((2*probs-1)*pgamma(lambda/e,alpha-1)+1-probs)/((2*probs-1)*pgamma(lambda/e,alpha)+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
