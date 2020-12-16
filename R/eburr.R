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
  
   fy=function(y){
     return((1+y^(-rho/gamma))^(1/rho))
   } 
   
   find_root=function(tau){
      fun=function(e){
        psi1=integrate(fy,e,Inf)$value
        return(psi1/(2*psi1+e+1/rho*beta(-1/rho+gamma/rho,1-gamma/rho))-1+tau)
      }
     
     return(uniroot(fun,c(0,100000))$root)
   }
   
   return(mapply(find_root,probs))
   
}