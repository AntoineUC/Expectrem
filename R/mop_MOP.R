mop.MOP=function (osx, k, p) 
{
  n = length(osx)
  dk = length(k)
  dp = length(p)
  km = max(k)
  nk = km + 1
  est = matrix(NA, dk, dp)
  for (j in 1:dp) if (p[j] == 0) {
    losx = rev(log(osx))
    losx = losx[1:nk]
    tmp = cumsum(losx)/(1:nk)
    estc = tmp[-nk] - losx[-1]
    est[, j] = estc[k]
  }
  else {
    tosx = rev(osx)
    tosx = (tosx[1:nk])^p[j]
    tmp = cumsum(tosx)/(1:nk)
    estc = tmp[-nk]/tosx[-1]
    estc = (1 - estc^-1)/p[j]
    est[, j] = estc[k]
  }
  return(list(EVI = est))
}