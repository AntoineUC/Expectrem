ESgumbel=function(probs, mu = 0, beta = 1){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(mu) > 1) {
    stop("mu must be of length 1.")
  }
  
  if (length(beta) > 1) {
    stop("beta must be of length 1.")
  }
  
  if(beta<0){
    stop("sigma must be positive.")
  }
  
  f=function(tau){
    -1/(1-tau)*sum(log(tau)^(1:100)/((1:100)*factorial(1:100)))-log(-log(tau))
  }
  
  return(mu+beta*sapply(probs,f))
}