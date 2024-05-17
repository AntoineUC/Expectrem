ESlogis=function(probs, mu = 0, sigma = 1){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(mu) > 1) {
    stop("mu must be of length 1.")
  }
  
  if (length(sigma) > 1) {
    stop("sigma must be of length 1.")
  }
  
  if(sigma<0){
    stop("sigma must be positive.")
  }
  
  return(mu-sigma*(probs/(1-probs)*log(probs)+log(1-probs)))
  
}