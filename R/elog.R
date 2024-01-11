elog=function(probs, mu = 0, sigma = 1, tol=1e-08,maxiter=100){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(sigma<0){
    stop("sigma must be positive.")
  }
  
  e=rep(0,length(probs))
  
  gap=1
  i=1
  while(gap >= tol && i<=maxiter){
    
    e1=(1-2*probs[which(probs*(1-probs) != 0)])*((exp(e[which(probs*(1-probs) != 0)])+1)*log(exp(e[which(probs*(1-probs) != 0)])+1)-e[which(probs*(1-probs) != 0)]*exp(e[which(probs*(1-probs) != 0)]))/((probs[which(probs*(1-probs) != 0)]-1)*exp(e[which(probs*(1-probs) != 0)])-probs[which(probs*(1-probs) != 0)])
    
    gap=max(abs(e1-e[which(probs*(1-probs) != 0)]),na.rm=T)
    
    e[which(probs*(1-probs) != 0)]=e1

    i=i+1
  }
  
  e[which(probs==0)]=-Inf
  
  e[which(probs==1)]=Inf

  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  
  return(mu+sigma*e)
  
}
