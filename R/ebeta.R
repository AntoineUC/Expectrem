ebeta=function (probs, shape1, shape2, niter = 50) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (shape1 <= 0 || shape2 <=0) {
    stop("shape parameters must be strictly positive.")
  }
  e = rep(shape1/(shape1+shape2), length(probs))
  i = 1
  while (i <= niter) {
    e = shape1/(shape1+shape2)*((2*probs-1)*(1-pbeta(e,shape1+1,shape2))+1-probs)/((2*probs-1)*(1-pbeta(e,shape1,shape2))+1-probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = 1
  return(e)
}