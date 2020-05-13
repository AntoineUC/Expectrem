tindexp=function(X,k,br=FALSE){

  if(!is.logical(br)){
    stop("br must be boolean.")
  }

  n=length(X)

  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")

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
