egpd=function(probs,gamma=0.5){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(gamma<=0 || gamma>=1){
    stop("gamma must be strictly between 0 and 1.")
  }
  
  fy=function(y){
    return((1+y*gamma)^(-1/gamma))
  } 
  
  find_root=function(tau){
    fun=function(e){
      psi1=integrate(fy,e,Inf)$value
      return(psi1/(2*psi1+e-1/(1-gamma))-1+tau)
    } 
    
    return(uniroot(fun,c(0,100000))$root)
  }
  
  return(mapply(find_root,probs))
  
}