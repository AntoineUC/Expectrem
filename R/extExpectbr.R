extExpectbr=function(X,k,tau,estim="Hill"){
  n=length(X)
  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")
  if(estim=="Hillbr"){
    gammahat=mopest$EVI[k]
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

  qtp=expect(X,1-k/n)

  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammahat*Fbar(X,qtp)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)

  rbet=(1-mean(X)/(qtp*(k/(n*(1-tau)))^(gammahat)))*(1/(2*tau-1))*(1+mopest$beta*gammahat*(1/gammahat-1)^(-mopest$rho)*(1-tau)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)

  return(qtp*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho)*(r/rbet)^gammahat*(1+((1/gammahat-1)^(-mopest$rho)*rbet^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(1-tau)^(-mopest$rho))/(1+((1/gammahat-1)^(-mopest$rho)*r^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(k/n)^(-mopest$rho)))
}
