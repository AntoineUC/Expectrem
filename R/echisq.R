echisq=function(probs,df,niter=50){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(df<=0){
    stop("df must be strictly greater than 0.")
  }
  
  e = rep(df, length(probs))
  i = 1
  while (i <= niter) {
    e = df*((2*probs-1)*(1-pgamma(e/2,df/2+1))+1-probs)/((2*probs-1)*(1-pgamma(e/2,df/2))+1-probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  
  return(e)
  
}