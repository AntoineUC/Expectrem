lpindex=function(X,k,p,br=FALSE){
  
  if(!is.logical(br)){
    stop("br must be boolean.")
  }
  
  n=length(X)
  
  if(any(k>n-1) || any(k<1)){
    stop("k must be between 1 and n-1.")
  }
  
  if(p<=1){
    stop("p must be greater than 1.")
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
    return(uniroot(func,c(0,10000000))$root)
  }
  
  if(br==FALSE){
    
    find_root=function(ks){
      
      ratio=mapply(Fbarhat1,mapply(quantp,1-ks/n))*n/ks
      
      f=function(gamma){
        return(gamma/beta(p,1/gamma-p+1))
      }
      
      func=function(y){
        return(f(y)-ratio)
      }
      return(uniroot(func,c(0.000001,1/(p-1)))$root)
    }
    
  }
  
  if(br==TRUE){
    
    mopest=mop(X,p=0,k=1:(n-1),method = "RBMOP")
    
    rhohat=mopest$rho
    
    bhat=mopest$beta
    
    meanlp=function(ks){
      return(mean(abs(X/quantp(1-ks/n)-1)^(p-1)))
    }
    
    find_root=function(ks){
      
      gammhat=mopest$EVI[ks]
      
      qp=mapply(quantp,1-ks/n)
      
      gpbar=gammhat/beta(p,(1/gammhat-p+1)*(1/gammhat-p+1>0))
      
      Kbar=gpbar^(-rhohat)/(gammhat^2*rhohat)*((1-rhohat)*beta(p,(1-rhohat)/gammhat-p+1)-beta(p,1/gammhat-p+1))
      
      matX=((matrix(rep(X,length(ks)),length(ks),n,byrow=T)-qp)<0)*((matrix(rep(X,length(ks)),length(ks),n,byrow=T)+qp)>0)
      
      ratio=mapply(Fbarhat1,qp)*n/ks*(1+bhat*gammhat*mapply(Fbarhat1,qp)^(-rhohat)*Kbar*gpbar^(1+rhohat))/mapply(meanlp,ks)
      
      ratio[which(ratio=='NaN')]=0
      
      f=function(gamma){
        return(gamma/beta(p,1/gamma-p+1))
      }
      
      func=function(y){
        return(f(y)-ratio)
      }
      return(uniroot(func,c(0.000001,1/(p-1)))$root)
    }
    
  }
  
  gammahat=mapply(find_root,k)
  
  return(gammahat)
}
