eweibull=function (probs, k, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL) 
{
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (k <= 0) {
    stop("k must be greater than 0.")
  }
  
  if (start.pt != "mean" && start.pt != "quantile" && start.pt != "custom") {
    stop("start.pt may be either custom, mean or quantile.")
  }
  
  if(start.pt=="mean"){
    e=rep(gamma(1 + 1/k), length(probs))
  }
  
  if(start.pt=="quantile"){
    e=qweibull(probs,shape=k)
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
    e1 = gamma(1 + 1/k) * ((2 * probs - 1) * (1-pgamma(e^k, 1 + 1/k)) + 1 - probs)/((2 * probs - 1) * exp(-e^k) + 
                                                                                      1 - probs)
    gap=max(abs(e1-e),na.rm=T)
    
    e=e1
    
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
