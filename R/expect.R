
expect=function (X, probs) {
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (is.vector(X)==FALSE) {
    stop("X must be a vector.")
  }
  
  if (sum(is.na(X))>0) {
    stop("X contains NA values.")
  }
  
  num1=(0:(n-2))*sort(X)[1:(n-1)]-c(0,cumsum(sort(X))[1:(n-2)])
  
  denum1=num1+sum(X)-cumsum(sort(X))[1:(n-1)]-(n-1:(n-1))*sort(X)[1:(n-1)]
  
  blow=num1/denum1
  
  f=function(tau){
    i=which((tau>=blow)==FALSE)[1]-1
    return((tau*(sum(X)-sum(sort(X)[1:i]))+(1-tau)*sum(sort(X)[1:i]))/(tau*n-(2*tau-1)*i))
  }
  
  e=mapply(f,probs)
  
  return(e)
  
}
