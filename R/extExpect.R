extExpect=function(X,k="kopt",tau,estim="Hill",method="direct",br=FALSE){

  n=length(X)

  mopest=mop(X, 1:(n-1), 0, method ="RBMOP")

  if(k=="kopt" && estim=="Hill"){
    k=trunc(((1-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho)))
  }

  if(k=="kopt" && estim=="tindexp"){
    k=trunc(((1-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho)))
    k=min(trunc(((1/mopest$EVI[k]-1)^(2*mopest$rho-1)*(1-mopest$EVI[k]-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2*abs(1-2*mopest$EVI[k])))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho))),trunc(n/2)-1)
  }

  if(k>n-1 || k<1){
    stop("k must be between 1 and n-1.")
  }

  if(estim!="Hill" && estim!="tindexp"){
    stop("estim may be either Hill or tindexp.")
  }

  if(method!="direct" && method!="indirect"){
    stop("method may be either direct or indirect.")
  }

  if(!is.logical(br)){
    stop("br must be boolean.")
  }

  if(estim=="Hill" && br==TRUE){
    gammahat=mopest$EVI[k]
  }
  if(estim=="Hill" && br==FALSE){
    gammahat=mop(X, 1:(n-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k,br)
  }

  qtp=expect(X,1-k/n)

  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammahat*Fbar(X,qtp)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)

  rbet=(1-mean(X)/(qtp*(k/(n*(1-tau)))^(gammahat)))*(1/(2*tau-1))*(1+mopest$beta*gammahat*(1/gammahat-1)^(-mopest$rho)*(1-tau)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)

  if(br==FALSE && method=="direct"){
    return(qtp*(k/(n*(1-tau)))^gammahat)
  }

  if(br==FALSE && method=="indirect"){
    return((1/gammahat-1)^(-gammahat)*quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat)
  }

  if(br==TRUE && method=="direct"){
    return(qtp*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho)*(r/rbet)^gammahat*(1+((1/gammahat-1)^(-mopest$rho)*rbet^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(1-tau)^(-mopest$rho))/(1+((1/gammahat-1)^(-mopest$rho)*r^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(k/n)^(-mopest$rho)))
  }

  if(br==TRUE && method=="indirect"){
    return((1/gammahat-1)^(-gammahat)*quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho)*(1+((1/gammahat-1)^(-mopest$rho)*rbet^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(1-tau)^(-mopest$rho))/(rbet^gammahat))
  }
}
