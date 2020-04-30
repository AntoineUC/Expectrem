extQuant=function(X,k,tau,estim="Hill"){
  n=length(X)
  if(estim=="Hillbr"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")$EVI[k]
  }
  if(estim=="Hill"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k)
  }
  if(estim=="tindexpbr"){
    gammahat=tindexpbr(X,k)
  }
  return(quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat)
}
