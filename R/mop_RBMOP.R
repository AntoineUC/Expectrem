mop.RBMOP=function (osx, k, p) 
{
  n = length(osx)
  losx = log(osx)
  dk = length(k)
  dp = length(p)
  est = matrix(NA, dk, dp)
  H = mop.MOP(osx, k, p)
  H = H$EVI
  rhoest = mop.rho(osx)
  betaest = mop.beta(losx, rhoest)
  for (j in 1:dp) est[, j] = H[, j] * (1 - ((betaest * (1 - 
                                                          p[j] * H[, j]))/(1 - rhoest - p[j] * H[, j])) * (n/k)^(rhoest))
  return(list(EVI = est, rho = rhoest, beta = betaest))
}