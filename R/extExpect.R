extExpect=function(X,k,tau,estim="Hill",br=FALSE){

  n=length(X)

  if(k>n-1 || k<1){
    stop("k must be between 1 and n-1.")
  }

  if(estim!="Hill" && estim!="tindexp"){
    stop("estim may be either Hill or tindexp.")
  }

  if(!is.logical(br)){
    stop("br must be boolean.")
  }

  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")

  if(estim=="Hill" && br==TRUE){
    gammahat=mopest$EVI[k]
  }
  if(estim=="Hill" && br==FALSE){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k,br)
  }

  qtp=expect(X,1-k/n)

  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammahat*Fbar(X,qtp)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)

  rbet=(1-mean(X)/(qtp*(k/(n*(1-tau)))^(gammahat)))*(1/(2*tau-1))*(1+mopest$beta*gammahat*(1/gammahat-1)^(-mopest$rho)*(1-tau)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)

  if(br==FALSE){
  return(qtp(k/(n*(1-tau)))^gammahat)
  }

  if(br==TRUE){
    return(qtp*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho)*(r/rbet)^gammahat*(1+((1/gammahat-1)^(-mopest$rho)*rbet^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(1-tau)^(-mopest$rho))/(1+((1/gammahat-1)^(-mopest$rho)*r^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(k/n)^(-mopest$rho)))
  }
}
