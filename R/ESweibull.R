ESweibull=function (probs, k) 
{
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(k) > 1) {
    stop("k must be of length 1.")
  }
  
  if (k <= 0) {
    stop("k must be greater than 0.")
  }
  
  return(gamma(1+1/k)*(1-pgamma(-log(1-probs),1+1/k))/(1-probs))
}