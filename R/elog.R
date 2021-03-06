elog=function(probs, mu = 0, sigma = 1, niter = 20){

  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }

  if(sigma<0){
    stop("sigma must be positive.")
  }

  e=rep(0,length(probs))

  i=1
  while(i<=niter){

    e=(1-2*probs)*((exp(e)+1)*log(exp(e)+1)-e*exp(e))/((probs-1)*exp(e)-probs)

    i=i+1
  }

  e[which(probs==0)]=-Inf

  e[which(probs==1)]=Inf

  return(mu+sigma*e)

}
