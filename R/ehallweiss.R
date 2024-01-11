ehallweiss=function (probs, alpha, beta, tol=1e-08, maxiter=100) {
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (alpha <= 1 || beta<=0) {
    stop("alpha and beta must be greater than 1 and 0, respectively.")
  }
  e = rep((2*(alpha-1)*(alpha+beta)+beta)/(2*(alpha-1)*(alpha+beta-1)), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<=maxiter) {
    e1 = ((2*probs-1)*(alpha/(alpha-1)*e^(1-alpha)+(alpha+beta)/(alpha+beta-1)*e^(1-alpha-beta))+(1-probs)*(2*(alpha-1)*(alpha+beta)+beta)/((alpha-1)*(alpha+beta-1)))/((2*probs-1)*(e^(-alpha)+e^(-alpha-beta))+2*(1-probs))
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  e[which(probs == 0)] = 1
  e[which(probs == 1)] = Inf
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
