ehallweiss=function (probs, alpha, beta, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL) {
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (alpha <= 1 || beta<=0) {
    stop("alpha and beta must be greater than 1 and 0, respectively.")
  }
  
  if (start.pt != "mean" && start.pt != "quantile" && start.pt != "custom") {
    stop("start.pt may be either custom, mean or quantile.")
  }
  
  if(start.pt=="mean"){
    e=rep((2*(alpha-1)*(alpha+beta)+beta)/(2*(alpha-1)*(alpha+beta-1)), length(probs))
  }
  
  if(start.pt=="quantile"){
    e=(2*(1-probs))^(-1/alpha)
  }
  
  if(start.pt=="custom"){
    if(length(x0) != length(probs)){
      stop("x0 and probs must have the same length.")
    }
    e=x0
  }
  
  gap=1
  i=1
  while (gap >= tol && i<=maxiter) {
    e1 = ((2*probs-1)*(alpha/(alpha-1)*e^(1-alpha)+(alpha+beta)/(alpha+beta-1)*e^(1-alpha-beta))+(1-probs)*(2*(alpha-1)*(alpha+beta)+beta)/((alpha-1)*(alpha+beta-1)))/((2*probs-1)*(e^(-alpha)+e^(-alpha-beta))+2*(1-probs))
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
