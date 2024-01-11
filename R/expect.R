
expect=function (X, probs, tol=1e-08, maxiter=100) {
 
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (is.vector(X)==FALSE) {
    stop("X must be a vector.")
  }
  
  if (sum(is.na(X))>0) {
    stop("X contains NA values.")
  }
   
  e = rep(mean(X), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<=maxiter) {
    e1 = ((2*probs-1)*apply(matrix(rep(X,length(probs)),length(X),length(probs))*((t(t(matrix(rep(X,length(probs)),length(X),length(probs)))-e))>0),2,mean)+(1-probs)*rep(mean(X),length(probs)))/((2*probs-1)*apply(((t(t(matrix(rep(X,length(probs)),length(X),length(probs)))-e))>0),2,mean)+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=i+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
  
}
