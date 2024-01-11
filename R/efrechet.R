efrechet=function(probs,a,tol=1e-08,maxiter=100){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(a<=1){
    stop("Expectiles do not exist if a is less than 1.")
  }
  
  e = rep(gamma(1-1/a), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<= maxiter) {
    e1 = gamma(1-1/a)*((2*probs[which(probs*(1-probs) != 0)]-1)*pgamma(e[which(probs*(1-probs) != 0)]^(-a),1-1/a)+1-probs[which(probs*(1-probs) != 0)])/((2*probs[which(probs*(1-probs) != 0)]-1)*(1-exp(-e[which(probs*(1-probs) != 0)]^(-a)))+1-probs[which(probs*(1-probs) != 0)])
    gap=max(abs(e1-e[which(probs*(1-probs) != 0)]),na.rm=T)
    e[which(probs*(1-probs) != 0)]=e1
    i=i+1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
