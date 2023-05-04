eburr=function(probs,a,b,niter=50){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(a<=0 || b<=0){
    stop("Parameters must be strictly positive.")
  }
  
  if(a*b<=1){
    stop("Expectiles do not exist if ab is less than 1.")
  }
  
  e = rep(b*beta(b-1/a,1/a+1), length(probs))
  i = 1
  while (i <= niter) {
    e = b*beta(b-1/a,1/a+1)*((2*probs-1)*pbeta(1/(1+e^a),b-1/a,1/a+1)+1-probs)/((2*probs-1)*(1+e^a)^(-b)+1-probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  
  return(e)
}
