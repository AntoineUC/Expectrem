egpd=function(probs,alpha,niter=50){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(alpha<=1){
    stop("alpha must be strictly greater than 1.")
  }
  
  e = rep(alpha/(alpha - 1), length(probs))
  i = 1
  while (i <= niter) {
    e = alpha/(alpha - 1) * ((2 * probs - 1) * e^(1 - alpha) + 
                               1 - probs)/((2 * probs - 1) * e^(-alpha) + 1 - probs)
    i = i + 1
  }
  e[which(probs == 0)] = 1
  e[which(probs == 1)] = Inf
  
  return((e-1)*alpha)
  
}
