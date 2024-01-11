et=function(probs, df, mu = 0, sigma = 1, tol=1e-08, maxiter=100){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(sigma<0){
    stop("sigma must be positive.")
  }
  
  if(df<=1){
    stop("df must greater than 1.")
  }
  
  e=rep(0,length(probs))
  
  gap=1
  i=1
  while(gap >= tol && i<=maxiter){
    
    e1=(2*probs-1)*df/(1-df)*dt(e,df)*(1+e^2/df)/((2*probs-1)*pt(e,df)-probs)
    
    gap=max(abs(e1-e),na.rm=T)
    
    e=e1

    i=i+1
  }

  if(i>maxiter){
    warning("Warning: maximum of iterations reached !")
  }
  
  return(mu+sigma*e)
}
