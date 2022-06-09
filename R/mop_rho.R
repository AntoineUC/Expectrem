mop.rho=function (osx) 
{
  losx = log(osx)
  n = length(losx)
  krho = floor(n^(0.995)):floor(n^(0.999))
  krhom = max(krho)
  nkrho = krhom + 1
  M = matrix(NA, length(krho), 3)
  tau0 = numeric(length(krho))
  tau1 = numeric(length(krho))
  W0 = numeric(length(krho))
  W1 = numeric(length(krho))
  losx = rev(losx)
  losx = losx[1:nkrho]
  estc1 = mop.MOP(osx, krho, p = 0)
  M[, 1] = estc1$EVI
  c12 = cumsum(losx^2)/(1:nkrho)
  c22 = (losx)^2
  c32 = cumsum(losx)/(1:nkrho)
  estc2 = c12[-nkrho] + c22[-1] - 2 * c32[-nkrho] * losx[-1]
  M[, 2] = estc2[krho]
  c13 = cumsum(losx^3)/(1:nkrho)
  c23 = (losx)^3
  estc3 = c13[-nkrho] - c23[-1] - 3 * (losx[-1] * c12[-nkrho] - 
                                         c32[-nkrho] * c22[-1])
  M[, 3] = estc3[krho]
  W0 = (log(M[, 1]) - (1/2) * log(M[, 2]/2))/((1/2) * log(M[, 
                                                            2]/2) - (1/3) * log(M[, 3]/6))
  W1 = (M[, 1] - (M[, 2]/2)^(1/2))/((M[, 2]/2)^(1/2) - (M[, 
                                                          3]/6)^(1/3))
  tau0 = -abs(3 * (W0 - 1)/(W0 - 3))
  tau1 = -abs(3 * (W1 - 1)/(W1 - 3))
  tau0median = median(tau0)
  tau1median = median(tau1)
  sumtau0 = as.numeric((tau0 - tau0median) %*% (tau0 - tau0median))
  sumtau1 = as.numeric((tau1 - tau1median) %*% (tau1 - tau1median))
  if ((sumtau0 < sumtau1) || (sumtau0 == sumtau1)) {
    return(tau0[length(krho)])
  }
  else {
    return(tau1[length(krho)])
  }
}