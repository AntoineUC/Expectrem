extQuant=function(X,k,tau,estim="Hill",br=FALSE){

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
  if(br==FALSE){
    return(quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat)
  }
  if(br==TRUE){
    return(quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho))
  }
}
