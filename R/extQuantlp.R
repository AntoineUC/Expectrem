extQuantlp=function(X,tau,k,p,estim="lpindex",br=FALSE){
  
  n=length(X)
  
  mopest=mop(X, 1:(n-1), 0, method ="RBMOP")
  
  if(any(k>n-1) || any(k<1)){
    stop("k must be between 1 and n-1.")
  }
  
    if(tau<=0 || tau>=1){
    stop("tau must be between 0 and 1.")
  }
  
  if(estim!="Hill" && estim!="lpindex"){
    stop("estim may be either Hill or lpindex.")
  }
  
  if(!is.logical(br)){
    stop("br must be boolean.")
  }
  
  Fbarhat1=function(y){
    return(sum(X>y)/n)
  }
  
  Fbarhatp=function(y){
    return(sum(abs(X-y)^(p-1)*(X>y))/sum(abs(X-y)^(p-1)))
  }
  
  quantp=function(alpha){
    func=function(y){
      return(Fbarhatp(y)-1+alpha)
    }
    return(uniroot(func,c(-100,10000000))$root)
  }
  
  meanl2=function(k){
    return(mean(abs(X/extExpect(X,tau=tau,k=k)-1)))
  }
  
  meanlp=function(k){
    return(mean(abs(X/mapply(quantp,1-k/n)-1)^(p-1)))
  }
  
  if(estim=="Hill" && br==TRUE){
    gammahat=mopest$EVI[k]
  }
  
  if(estim=="Hill" && br==FALSE){
    gammahat=mop(X, 1:(n-1), 0, method ="MOP")$EVI[k]
  }
  
  if(estim=="lpindex"){
    gammahat=lpindex(X,k,p,br)
  }
  
  gphat=gammahat/beta(p,1/gammahat-p+1)
  
  if(br==FALSE){
    
    quantilp=((1-tau)*n/k)^(-gammahat)*mapply(quantp,1-k/n)*(gphat)^gammahat
    
  }
  
  if(br==TRUE){
    
    rhohat=mopest$rho
    
    bhat=mopest$beta
    
    B1=1+(((1-tau)*n/k)^(-rhohat)-1)/rhohat*bhat*gammahat*(n/k)^(rhohat)
    
    Kphat=gphat^(-rhohat)/(gammahat^2*rhohat)*((1-rhohat)*beta(p,(1-rhohat)/gammahat-p+1)-beta(p,1/gammahat-p+1))
    
    rp=(1+bhat*gammahat*mapply(Fbarhat1,mapply(quantp,1-k/n))^(-rhohat)*Kphat*gphat^(1+rhohat))^(-1)*mapply(meanlp,k)
    
    B3=1+(gphat^(-rhohat)*rp^(-rhohat)-1)/rhohat*bhat*gammahat*(k/n)^(-rhohat)
    
    quantilp=((1-tau)*n/k)^(-gammahat)*mapply(quantp,1-k/n)*(gphat)^gammahat*B1*1/B3*(rp)^gammahat
    
  }
  
  return(quantilp)
}
