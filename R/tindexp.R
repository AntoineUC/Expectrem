tindexp=function(X,k="kopt",br=FALSE){
  
  if(!is.logical(br)){
    stop("br must be boolean.")
  }
  
  n=length(X)
  
  if(any(k>n-1) || any(k<1)){
    stop("k must be between 1 and n-1.")
  }
  
  mopest=mop(X, 1:(length(X)-1), 0, method ="RBMOP")
  
  if(k[1]=="kopt"){
    k=trunc(((1-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho)))
    k=min(trunc(((1/mopest$EVI[k]-1)^(2*mopest$rho-1)*(1-mopest$EVI[k]-mopest$rho)^2/(-2*mopest$rho*mopest$beta^2*abs(1-2*mopest$EVI[k])))^(1/(1-2*mopest$rho))*n^(-2*mopest$rho/(1-2*mopest$rho))),trunc(n/2)-1)
  }
  
  #gammhat=mopest$EVI[k] #previous version, published in Girard, S., Stupfler, G. and Usseglio-Carleve, A., On automatic bias reduction for extreme expectile estimation, Statistics and Computing, 32(64).
  
  qtp=expect(X,1-k/n)
  
  gammhat = 1/(1 + Fbar(X, qtp)/(k/n)) #new version (2023)
  
  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammhat*Fbar(X,qtp)^(-mopest$rho)/(gammhat*(1-mopest$rho-gammhat)))^(-1)
  
  if(br==FALSE){
    return(1/(1+Fbar(X,qtp)/(k/n)))
  }
  if(br==TRUE){
    return(1/(1+Fbar(X,qtp)*n/k*1/r))
  }
}
