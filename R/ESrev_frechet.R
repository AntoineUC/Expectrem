ESrev_frechet=function (probs, a) 
{
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(a) > 1) {
    stop("k must be of length 1.")
  }
  
  if (a <= 0) {
    stop("a must be greater than 0.")
  }
  
  return(1-(pgamma(-log(probs),1/a)*gamma(1/a)-a*probs*(-log(probs))^(1/a))/(a*(1-probs)))
}