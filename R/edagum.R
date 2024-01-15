edagum=function(probs, a, b, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(b<=0){
    stop("b must be strictly positive.")
  }
  
  if(a<=1){
    stop("Expectiles do not exist if a is less than 1.")
  }
  
  if (start.pt != "mean" && start.pt != "quantile" && start.pt != "custom") {
    stop("start.pt may be either custom, mean or quantile.")
  }
  
  if(start.pt=="mean"){
    e=rep(b*beta(b+1/a,1-1/a), length(probs))
  }
  
  if(start.pt=="quantile"){
    e=(probs^(-1/b)-1)^(-1/a)
  }
  
  if(start.pt=="custom"){
    if(length(x0) != length(probs)){
      stop("x0 and probs must have the same length.")
    }
    e=x0
  }
  
  gap=1
  i=1
  while (gap >= tol && i<= maxiter) {
    e1 = b*beta(b+1/a,1-1/a)*((2*probs-1)*(1-pbeta(1/(1+e^(-a)),b+1/a,1-1/a))+1-probs)/(probs-(2*probs-1)*(1+e^(-a))^(-b))
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
