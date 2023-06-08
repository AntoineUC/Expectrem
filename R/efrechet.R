efrechet=function(probs,a,niter=50){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(a<=1){
    stop("Expectiles do not exist if a is less than 1.")
  }
  
  e = rep(gamma(1-1/a), length(probs))
  i = 1
  while (i <= niter) {
    e = gamma(1-1/a)*((2*probs-1)*pgamma(e^(-a),1-1/a)+1-probs)/((2*probs-1)*(1-exp(-e^(-a)))+1-probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  
  return(e)
}