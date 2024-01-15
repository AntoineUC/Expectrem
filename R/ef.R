ef=function (probs, df1, df2, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL) 
{
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (df2 <= 2) {
    stop("df2 must be strictly greater than 2.")
  }
  
  if (df1 <= 0) {
    stop("df1 must be strictly greater than 0.")
  }
  
  if (start.pt != "mean" && start.pt != "quantile" && start.pt != "custom") {
    stop("start.pt may be either custom, mean or quantile.")
  }
  
  if(start.pt=="mean"){
    e=rep(df2/(df2-2), length(probs))
  }
  
  if(start.pt=="quantile"){
    e=qf(probs,df1,df2)
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
    e1 = df2/(df2-2)*((2*probs-1)*(1-pbeta(e*df1/(df2+e*df1),df1/2+1,df2/2-1))+1-probs)/((2*probs-1)*(1-pbeta(e*df1/(df2+df1*e),df1/2,df2/2))+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=+1
  }
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
