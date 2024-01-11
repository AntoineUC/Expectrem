ef=function (probs, df1, df2, tol=1e-08,maxiter=100) 
{
  if (min(probs) < 0 || max(probs) > 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (df2 <= 2) {
    stop("df2 must be strictly greater than 2.")
  }
  if (df1 <= 0) {
    stop("df1 must be strictly greater than 0.")
  }
  e = rep(df2/(df2-2), length(probs))
  gap=1
  i=1
  while (gap >= tol && i<= maxiter) {
    e1 = df2/(df2-2)*((2*probs-1)*(1-pbeta(e*df1/(df2+e*df1),df1/2+1,df2/2-1))+1-probs)/((2*probs-1)*(1-pbeta(e*df1/(df2+df1*e),df1/2,df2/2))+1-probs)
    gap=max(abs(e1-e),na.rm=T)
    e=e1
    i=+1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  return(e)
}
