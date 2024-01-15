elogis=function(probs, mu = 0, sigma = 1, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(sigma<0){
    stop("sigma must be positive.")
  }
  
  if (start.pt != "mean" && start.pt != "quantile" && start.pt != "custom") {
    stop("start.pt may be either custom, mean or quantile.")
  }
  
  if(start.pt=="mean"){
    e=rep(0,length(probs))
  }
  
  if(start.pt=="quantile"){
    e=qlogis(probs)
  }
  
  if(start.pt=="custom"){
    if(length(x0) != length(probs)){
      stop("x0 and probs must have the same length.")
    }
    e=x0
  }
  
  gap=1
  i=1
  while(gap >= tol && i<=maxiter){
    
    e1=(1-2*probs)*((exp(e)+1)*log(exp(e)+1)-e*exp(e))/((probs-1)*exp(e)-probs)
    
    gap=max(abs(e1-e),na.rm=T)
    
    e=e1
    
    i=i+1
  }
  
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  
  return(mu+sigma*e)
  
}
