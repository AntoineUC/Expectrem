ESkumaraswamy=function (probs, a, b) 
{
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (length(a) > 1) {
    stop("a must be of length 1.")
  }
  
  if (length(b) > 1) {
    stop("b must be of length 1.")
  }
  if (a <= 0 || b <=0) {
    stop("shape parameters must be strictly positive.")
  }
  
  return(b*beta(1/a+1,b)*(1-pbeta(1-(1-probs)^(1/b),1/a+1,b))/(1-probs))
}