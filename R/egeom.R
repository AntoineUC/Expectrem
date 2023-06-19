egeom=function (probs, p) 
{
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  if (p < 0 || p > 1) {
    stop("p must be between 0 and 1.")
  }
  e = rep(1/p, length(probs))
  xp=1/p-1/log(1-p)*lambertWp(-(1-p)^(1/p)*log(1-p)/p*(2*probs-1)/(1-probs))
  e=((2*probs-1)*(1-p)^(trunc(xp))*(1+p*trunc(xp))+1-probs)/(p*((2*probs-1)*(1-p)^(trunc(xp))+1-probs))-1
  return(e)
}
