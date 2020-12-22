eburr=function(probs,gamma=0.5,rho=-1){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(gamma<=0 || gamma>=1){
    stop("gamma must be strictly between 0 and 1.")
  }
  
  if(rho>0){
    stop("rho must be negative.")
  }
  
  kmax=100
  
  coeffs=cumprod((1/rho-0:(kmax-1)))/factorial(1:kmax)
  
  find_root=function(tau){
    fun=function(e){
      if(e>=1){
        psi1=-e^(1-1/gamma)/(1-1/gamma)-sum(coeffs*e^(1+(1:kmax)*rho/gamma-1/gamma)/(1+(1:kmax)*rho/gamma-1/gamma))
      }
      if(e<1){
        psi1=1-e+sum(coeffs*(1-e^(1-(1:kmax)*rho/gamma))/(1-(1:kmax)*rho/gamma))-1/(1-1/gamma)-sum(coeffs/(1+(1:kmax)*rho/gamma-1/gamma))
        
      }
      return(psi1/(2*psi1+e+1/rho*beta(-1/rho+gamma/rho,1-gamma/rho))-1+tau)
    }
    
    return(uniroot(fun,c(0,10000000))$root)
  }
  
  return(mapply(find_root,probs))
  
}
