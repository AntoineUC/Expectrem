echisq=function(probs,df,tol=1e-08,maxiter=100){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(df<=0){
    stop("df must be strictly greater than 0.")
  }
  
  e = rep(df, length(probs))
  gap=1
  i=1
  while (gap >= tol && i<= maxiter) {
    e1 = df*((2*probs-1)*(1-pgamma(e/2,df/2+1))+1-probs)/((2*probs-1)*(1-pgamma(e/2,df/2))+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
  
}
