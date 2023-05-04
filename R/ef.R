ef=function (probs, df1, df2, niter = 50) 
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
  i = 1
  while (i <= niter) {
    e = df2/(df2-2)*((2*probs-1)*(1-pbeta(e*df1/(df2+e*df1),df1/2+1,df2/2-1))+1-probs)/((2*probs-1)*(1-pbeta(e*df1/(df2+df1*e),df1/2,df2/2))+1-probs)
    i = i + 1
  }
  e[which(probs == 0)] = 0
  e[which(probs == 1)] = Inf
  return(e)
}