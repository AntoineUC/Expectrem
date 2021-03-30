extExpectlp=function(X,tau,k,p,estim="lpindex",br=FALSE){
  
  n=length(X)
  
  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")
  
  if(k>n-1 || k<1){
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
    return(uniroot(func,c(-100,100000))$root)
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
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  
  if(estim=="lpindex"){
    gammahat=lpindex(X[which(X>0)],k,p,br)
  }
  
  g2hat=1/gammahat-1
  
  gphat=gammahat/beta(p,1/gammahat-p+1)
  
  if(br==FALSE){
    
    expectilp=((1-tau)*n/k)^(-gammahat)*mapply(quantp,1-k/n)*(gphat/g2hat)^gammahat
    
  }
  
  if(br==TRUE){
    
    rhohat=mopest$rho
    
    bhat=mopest$beta
    
    B1=1+(((1-tau)*n/k)^(-rhohat)-1)/rhohat*bhat*gammahat*(n/k)^(rhohat)
    
    K2hat=g2hat^(-rhohat)/(gammahat^2*rhohat)*((1-rhohat)*beta(2,(1-rhohat)/gammahat-1)-beta(2,1/gammahat-1))
    
    Kphat=gphat^(-rhohat)/(gammahat^2*rhohat)*((1-rhohat)*beta(p,(1-rhohat)/gammahat-p+1)-beta(p,1/gammahat-p+1))
    
    r2=(1+bhat*gammahat*mapply(Fbarhat1,extExpect(X,tau=tau,k=k))^(-rhohat)*K2hat*g2hat^(1+rhohat))^(-1)*mapply(meanl2,k)
    
    rp=(1+bhat*gammahat*mapply(Fbarhat1,mapply(quantp,1-k/n))^(-rhohat)*Kphat*gphat^(1+rhohat))^(-1)*mapply(meanlp,k)
    
    B2=1+(g2hat^(-rhohat)*r2^(-rhohat)-1)/rhohat*bhat*gammahat*(1-tau)^(-rhohat)
    
    B3=1+(gphat^(-rhohat)*rp^(-rhohat)-1)/rhohat*bhat*gammahat*(k/n)^(-rhohat)
    
    expectilp=((1-tau)*n/k)^(-gammahat)*mapply(quantp,1-k/n)*(gphat/g2hat)^gammahat*B1*B2/B3*(rp/r2)^gammahat
    
  }
  
  return(expectilp)
}
