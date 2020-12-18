epareto=function(probs,gamma=0.5){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(gamma<=0 || gamma>=1){
    stop("gamma must be strictly between 0 and 1.")
  }
  
  find_root=function(tau){
    fun=function(e){
      psi1=-e^(1-1/gamma)/(1-1/gamma)
      return(psi1/(2*psi1+e-1/(1-gamma))-1+tau)
    } 
    
    return(uniroot(fun,c(1,10000000))$root)
  }
  
  mapply(find_root,probs)
  

}