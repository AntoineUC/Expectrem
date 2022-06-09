mop.beta=function (losx, rhoest) 
{
  n = length(losx)
  k1 = floor(n^(0.999))
  v = c(0, rhoest, 2 * rhoest)
  D = numeric(length(v))
  c1 = ((1:k1)/k1)^(-v[1])
  c2 = ((1:k1)/k1)^(-v[2])
  c3 = ((1:k1)/k1)^(-v[3])
  c = (1:k1) * ((losx[n:(n - k1 + 1)]) - (losx[(n - 1):(n - 
                                                          k1)]))
  D[1] = sum(c1 * c)/k1
  D[2] = sum(c2 * c)/k1
  D[3] = sum(c3 * c)/k1
  d = sum(c2)/k1
  betaest = (k1/n)^(v[2]) * (d * D[1] - D[2])/(d * D[2] - D[3])
  return(betaest)
}