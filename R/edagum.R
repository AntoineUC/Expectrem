edagum=function(probs,a,b,niter=50){
  
  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(b<=0){
    stop("b must be strictly positive.")
  }
  
  if(a<=1){
    stop("Expectiles do not exist if a is less than 1.")
  }
  
  e = rep(b*beta(b+1/a,1-1/a), length(probs))
  i = 1
  while (i <= niter) {
    e = b*beta(b+1/a,1-1/a)*((2*probs-1)*(1-pbeta(1/(1+e^(-a)),b+1/a,1-1/a))+1-probs)/(probs-(2*probs-1)*(1+e^(-a))^(-b))
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  
}