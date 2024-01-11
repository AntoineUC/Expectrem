enorm=function(probs, mu = 0, sigma = 1, tol=1e-08, maxiter=100){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(sigma<0){
    stop("sigma must be positive.")
  }
  
  e=rep(0,length(probs))
  
  gap=1
  i=1
  while(gap >= tol && i<=maxiter){
    
    e1=(1-2*probs)*dnorm(e)/((2*probs-1)*pnorm(e)-probs)
    
    gap=max(abs(e1-e),na.rm=T)
    
    e=e1

    i=i+1
  }

  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  
  return(mu+sigma*e)
}
