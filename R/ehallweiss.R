ehallweiss=function (probs, alpha=2, beta=1, niter = 50) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (alpha <= 1 || beta<=0) {
    stop("alpha and beta must be greater than 1 and 0, respectively.")
  }
  e = rep((2*(alpha-1)*(alpha+beta)+beta)/(2*(alpha-1)*(alpha+beta-1)), length(probs))
  i = 1
  while (i <= niter) {
    e = ((2*probs-1)*(alpha/(alpha-1)*e^(1-alpha)+(alpha+beta)/(alpha+beta-1)*e^(1-alpha-beta))+(1-probs)*(2*(alpha-1)*(alpha+beta)+beta)/((alpha-1)*(alpha+beta-1)))/((2*probs-1)*(e^(-alpha)+e^(-alpha-beta))+2*(1-probs))
    i = i + 1
  }
  e[which(probs == 0)] = 1
  e[which(probs == 1)] = Inf
  return(e)
}