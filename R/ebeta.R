ebeta=function (probs, shape1, shape2, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL) 
{
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (shape1 <= 0 || shape2 <=0) {
    stop("shape parameters must be strictly positive.")
  }
  if (start.pt != "mean" && start.pt != "quantile" && start.pt != "custom") {
    stop("start.pt may be either custom, mean or quantile.")
  }
  
  if(start.pt=="mean"){
    e=rep(shape1/(shape1 + shape2), length(probs))
  }
  
  if(start.pt=="quantile"){
    e=qbeta(probs,shape1,shape2)
  }
  
  if(start.pt=="custom"){
    if(length(x0) != length(probs)){
      stop("x0 and probs must have the same length.")
    }
    e=x0
  }

  gap=1
  i=1
  while (gap >= tol && i<= maxiter) {
    e1 = shape1/(shape1+shape2)*((2*probs-1)*(1-pbeta(e,shape1+1,shape2))+1-probs)/((2*probs-1)*(1-pbeta(e,shape1,shape2))+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
