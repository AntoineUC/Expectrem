et=function(probs, df = 2, mu = 0, sigma = 1, niter = 20){

  if (min(probs) < 0 || max(probs) > 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }

  if(sigma<0){
    stop("sigma must be positive.")
  }

  if(df<=1){
    stop("df must greater than 1.")
  }

  e=rep(0,length(probs))

  i=1
  while(i<=niter){

    e=(2*probs-1)*df/(1-df)*dt(e,df)*(1+e^2/df)/((2*probs-1)*pt(e,df)-probs)

    i=i+1
  }

  e[which(probs==0)]=-Inf

  e[which(probs==1)]=Inf

  return(mu+sigma*e)
}
