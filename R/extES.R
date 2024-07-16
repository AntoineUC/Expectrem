extES=function(X,k,tau,estim="Moment",method="direct",ci.level=0.95,n.bootstrap=10000){
  
  if (length(ci.level) > 1) {
    stop("ci.level must be of length 1.")
  }
  
  if (length(k) > 1) {
    stop("k must be of length 1.")
  }
  
  if(k>n-1 || k<1){
    stop("k must be between 1 and n-1.")
  }
  
  if (tau >= 1 || tau <= 0) {
    stop("tau must be strictly between 0 and 1.")
  }
  
  if(n.bootstrap<1){
    stop("n.bootstrap must be a positive integer.")
  }
  
  if (ci.level >= 1 || ci.level <= 0) {
    stop("ci.level must be strictly between 0 and 1.")
  }
  
  if(method!="direct" && method!="indirect"){
    stop("method may be either direct or indirect.")
  }
  
  if(estim!="Hill" && estim!="Moment"){
    stop("estim may be either Hill or Moment.")
  }
  
  n=length(X)
  
  quant=quantile(X,1-k/n,type=1)
  
  M1=mean(log(quantile(X,1-(0:(k-1))/n,type=1)/quantile(X,1-k/n,type=1)))
  
  M2=mean(log(quantile(X,1-(0:(k-1))/n,type=1)/quantile(X,1-k/n,type=1))^2)
  
  M=M1+1-0.5*(1-M1^2/M2)^(-1)
  
  ahat=quant*M1*(1-M+M1)
  
  Int1=((k/(n*(1-tau)))^M-1)/M
  
  if(method=="direct" && estim=="Moment"){
    
    set.seed(1234)
    
    Y=matrix((1-runif(n*n.bootstrap))^(-1),n,n.bootstrap)
    
    logY=log(Y)[1:k,]
    
    Y_kn=exp(apply(log(Y[((k+1):n),])/matrix(rep(((k+1):n),),n-k,n.bootstrap,byrow=F),2,sum))
    
    Z=runif(n.bootstrap)
    
    S2=exp(apply(log(Y[((k+1):n),])/matrix(rep(((k+1):n),),n-k,n.bootstrap,byrow=F),2,sum))*k/n
    
    varD=(M^2+1)*(M>=0)+((1-M)^2*(1-2*M)*(1-M+6*M^2))/((1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
    
    vara=(M^2+2)*(M>=0)+(2-16*M+51*M^2-69*M^3+50*M^4-24*M^5)/((1-2*M*(M<0))*(1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
    
    covgamma=(M-1)*(M>=0)+((1-M)^2*(-1+4*M-12*M^2))/((1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
    
    Sigma=matrix(c(vara,covgamma,covgamma,varD),2,2)
    
    gammahatbis=M+sqrt((varD)/k)*qnorm(Z*pnorm((1-M)/sqrt((varD)/k)))
    
    matgammahatbis=matrix(rep(gammahatbis,k),k,n.bootstrap,byrow=T)
    
    cteint=mean(X[which(X>quant)])
    
    I1=((k/(n*(1-tau)))^gammahatbis-1)/(gammahatbis)
    
    I2=((k/(n*(1-tau)))^gammahatbis*log(k/(n*(1-tau)))*gammahatbis-(k/(n*(1-tau)))^gammahatbis+1)/gammahatbis^2
    
    Y_p_mat=Y[1:k,]^matgammahatbis
    
    Y_k_vec=(Y_kn*k/n)^gammahatbis
    
    Gn1=1/(gammahatbis*abs(1-gammahatbis))-1/gammahatbis*apply(Y_p_mat,2,mean)*Y_k_vec #abs sur le 1-gamma
    
    Wbis=(Y_p_mat-1)/matgammahatbis*(matgammahatbis<0)+logY*(matgammahatbis>=0)
    
    rac1=apply(Wbis,2,mean)
    
    rac2=apply(Wbis^2,2,mean)
    
    Gn2=0.5/(1-rac1^2/rac2)*(rac1)+gammahatbis*(gammahatbis>0)*(Y_k_vec-1)/gammahatbis
    
    Gn3=gammahatbis*(gammahatbis>0)*(rac1-1)+(1-0.5/(1-rac1^2/rac2)-gammahatbis*(gammahatbis<0))
    
    Ln=Gn1/Gn2+I1/(1-gammahatbis)*(1/Gn2-1)-(I1/(1-gammahatbis)^2+I2/(1-gammahatbis))*Gn3
    
    ctehat=cteint+ahat/(1-M)*(Int1)
    
    Lup=quantile(Ln,1-(1-ci.level)/2,type=8)
    
    Ldown=quantile(Ln,(1-ci.level)/2,type=8)
    
    cteup=ctehat+ahat*Lup
    
    ctedown=ctehat+ahat*Ldown
    
  }
  
  if(method=="indirect" && estim=="Moment"){
    
    Int2=((k/(n*(1-tau)))^M*log(k/(n*(1-tau)))*M-(k/(n*(1-tau)))^M+1)/M^2
    
    Int3=((k/(n*(1-tau)))^M*log(k/(n*(1-tau)))^2*M^2-2*(k/(n*(1-tau)))^M*log(k/(n*(1-tau)))*M+2*(k/(n*(1-tau)))^M-2)/M^3
    
    ctehat=quant+ahat/(1-M)*(1+Int1)
    
    vgamma=(M^2+1)*(M>=0)+((1-M)^2*(1-2*M)*(1-M+6*M^2))/((1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
    
    va=(M^2+2)*(M>=0)+(2-16*M+51*M^2-69*M^3+50*M^4-24*M^5)/((1-2*M*(M<0))*(1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
    
    covagamma=(M-1)*(M>=0)+((1-M)^2*(-1+4*M-12*M^2))/((1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
    
    vaprime=2*M*(M>=0)+2*(1-M-53*M^2+337*M^3-1020*M^4+1732*M^5-1536*M^6+576*M^7)/((1-4*M)^2*(1-3*M)^2*(1-2*M)^2)*(M<0)
    
    vgammaprime=2*M*(M>=0)+2*(1+3*M-69*M^2+281*M^3-552*M^4+552*M^5-216*M^6)/(1-7*M+12*M^2)^2*(M<0)
    
    covagammaprime=(M>=0)+(-1-18*M+159*M^2-440*M^3+588*M^4-288*M^5)/(1-7*M+12*M^2)^2*(M<0)
    
    Lambda22=sqrt(va-M^2)
    
    Lambda23=covagamma/Lambda22
    
    Lambda33=sqrt(vgamma-covagamma^2/(va-M^2))
    
    Lambda=matrix(c(1,M,0,0,Lambda22,Lambda23,0,0,Lambda33),3,3)
    
    Lambda22prime=(vaprime-2*M)/(2*Lambda22)
    
    Lambda23prime=(covagammaprime*Lambda22-covagamma*(vaprime-2*M)/(2*Lambda22))/(va-M^2)
    
    Lambda33prime=(vgammaprime-(2*covagammaprime*covagamma*(va-M^2)-covagamma^2*(vaprime-2*M))/(va-M^2)^2)/(2*Lambda33)
    
    b0=c((1+M*Int1)/(1-M),(1+Int1)/(1-M)*Lambda22+Lambda23*(1+Int1+(1-M)*Int2)/(1-M)^2,Lambda33*(1+Int1+(1-M)*Int2)/(1-M)^2)
    
    a2=(1+Int1)/(1-M)
    
    a3=(1+Int1+(1-M)*Int2)/(1-M)^2
    
    a2prime=(1+Int1+(1-M)*Int2)/(1-M)^2
    
    a3prime=(2*(1+Int1+(1-M)*Int2)+(1-M)^2*Int3)/(1-M)^3
    
    theta1=a2+M*a2prime
    
    theta2=a2prime*Lambda22+a3prime*Lambda23+a2*Lambda22prime+a3*Lambda23prime
    
    theta3=a3prime*Lambda33+a3*Lambda33prime
    
    theta=c(theta1,theta2,theta3)
    
    A=matrix(c(0,-1/(2*sqrt(k)),0,-1/(2*sqrt(k)),-(1+Int1)/((1-M)*sqrt(k)),-(1+Int1+(1-M)*Int2)/(2*(1-M)^2*sqrt(k)),0,-(1+Int1+(1-M)*Int2)/(2*(1-M)^2*sqrt(k)),0),3,3)
    
    B=t(Lambda)%*%A%*%Lambda-matrix(c(0,theta1*Lambda23/(2*sqrt(k)),theta1*Lambda33/(2*sqrt(k)),theta1*Lambda23/(2*sqrt(k)),theta2*Lambda23/(sqrt(k)),(theta2*Lambda33+theta3*Lambda23)/(2*sqrt(k)),theta1*Lambda33/(2*sqrt(k)),(theta2*Lambda33+theta3*Lambda23)/(2*sqrt(k)),theta3*Lambda33/(sqrt(k))),3,3)
    
    omega=eigen(B)$values
    
    lambda=(t(eigen(B)$vectors)%*%b0)[,1]^2/(4*omega^2)
    
    m=-sum(omega*lambda)
    
    variance=2*sum(diag(B%*%B))+sum(b0^2)
    
    quantup=sqrt(variance)*qnorm(1-(1-ci.level)/2)-sum(diag(B))
    
    quantdown=-sqrt(variance)*qnorm(1-(1-ci.level)/2)-sum(diag(B))
    
    cteup=ctehat+ahat/sqrt(k)*quantup
    
    ctedown=ctehat+ahat/sqrt(k)*quantdown
    
  }
  
  if(method=="direct" && estim=="Hill"){
    
    cteint=mean(X[which(X>quant)])
    
    ctehat=(k/(n*(1-tau)))*M1*cteint
    
    ctedown=ctehat*exp(-M1*log(k/(n*(1-tau)))/sqrt(k)*qnorm(1-(1-ci.level)/2))
    
    cteup=ctehat*exp(M1*log(k/(n*(1-tau)))/sqrt(k)*qnorm(1-(1-ci.level)/2))
    
  }
  
  if(method=="indirect" && estim=="Hill"){
    
    ctehat=(k/(n*(1-tau)))*M1*quant/(1-M1)
    
    ctedown=ctehat*exp(-M1*log(k/(n*(1-tau)))/sqrt(k)*qnorm(1-(1-ci.level)/2))
    
    cteup=ctehat*exp(M1*log(k/(n*(1-tau)))/sqrt(k)*qnorm(1-(1-ci.level)/2))
    
  }
  
  return(list(Lower_bound=ctedown,Point_estimate=ctehat,Upper_bound=cteup))  
  
}