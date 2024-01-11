edagum=function(probs,a,b,tol=1e-08,maxiter=100){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(b<=0){
    stop("b must be strictly positive.")
  }
  
  if(a<=1){
    stop("Expectiles do not exist if a is less than 1.")
  }
  
  e = rep(b*beta(b+1/a,1-1/a), length(probs))
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
