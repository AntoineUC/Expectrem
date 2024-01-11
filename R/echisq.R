echisq=function(probs,df,tol=1e-08){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(df<=0){
    stop("df must be strictly greater than 0.")
  }
  
  e = rep(df, length(probs))
  gap=1
  while (gap >= tol) {
    e1 = df*((2*probs-1)*(1-pgamma(e/2,df/2+1))+1-probs)/((2*probs-1)*(1-pgamma(e/2,df/2))+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  
  return(e)
  
}
