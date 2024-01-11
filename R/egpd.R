egpd=function(probs,alpha,tol=1e-08){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(alpha<=1){
    stop("alpha must be strictly greater than 1.")
  }
  
  e = rep(alpha/(alpha - 1), length(probs))
  gap=1
  while (gap >= tol) {
    e1 = alpha/(alpha - 1) * ((2 * probs - 1) * e^(1 - alpha) + 
                               1 - probs)/((2 * probs - 1) * e^(-alpha) + 1 - probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
  }
  e[which(probs == 0)] = 1
  e[which(probs == 1)] = Inf
  
  return((e-1)*alpha)
  
}
