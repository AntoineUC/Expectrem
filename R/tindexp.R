tindexp=function(X,k="kopt",br=FALSE){

  if(!is.logical(br)){
    stop("br must be boolean.")
  }

  n=length(X)

  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")

  if(k[1]=="kopt"){
    k=trunc(((1-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho)))
    k=min(trunc(((1/mopest$EVI[k]-1)^(2*mopest$rho-1)*(1-mopest$EVI[k]-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2*abs(1-2*mopest$EVI[k])))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho))),trunc(n/2)-1)
  }

  gammhill=mopest$EVI[k]

  qtp=expect(X,1-k/n)

  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammhill*Fbar(X,qtp)^(-mopest$rho)/(gammhill*(1-mopest$rho-gammhill)))^(-1)

  if(br==FALSE){
    return(1/(1+Fbar(X,qtp)/(k/n)))
  }
  if(br==TRUE){
    return(1/(1+Fbar(X,qtp)*n/k*1/r))
  }
}
