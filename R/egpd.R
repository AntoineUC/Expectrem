egpd=function(probs,gamma=0.5,method="uniroot"){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(gamma<=0 || gamma>=1){
    stop("gamma must be strictly between 0 and 1.")
  }
  
  if(method!="mc" && method!="uniroot"){
    stop("method must be either mc or uniroot.")
  }
  
  find_root=function(tau){
    fun=function(e){
      psi1=(1+gamma*e)^(1-1/gamma)/(1-gamma)
      return(psi1/(2*psi1+e-1/(1-gamma))-1+tau)
    } 
    
    return(uniroot(fun,c(0,10000000))$root)
  }
  
  if(method=="uniroot"){
    return(mapply(find_root,probs))
  }
  
  if(method=="mc"){
    return(expect(((1-runif(10000000))^(-gamma)-1)/gamma,probs))
  }
}