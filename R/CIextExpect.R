CIextExpect=function(X, k = trunc(length(X)/10), tau, method = "direct",ci.level=0.95) 
{
  n = length(X)
  mopest = mop(X, 1:(n - 1), 0, method = "RBMOP")
  if (length(k) > 1) {
    stop("k must be of length 1.")
  }
  if (length(tau) > 1) {
    stop("tau must be of length 1.")
  }
  if (k > n - 1 || k < 1) {
    stop("k must be between 1 and n-1.")
  }
  if (method != "direct" && method != "indirect") {
    stop("method may be either direct or indirect.")
  }
  
  if (method == "direct") {
    qtp = expect(X, 1 - k/n)
    gammahat=tindexp(X,k,br=T)
    if(gammahat>0.5){
      print("WARNING : Tail index above 1/2 ! Use the indirect method rather ?")
    }
    r = (1 - mean(X)/qtp) * (n/(n - 2 * k)) * (1 + mopest$beta * 
                                                 gammahat * Fbar(X, qtp)^(-mopest$rho)/(gammahat * (1 - 
                                                                                                      mopest$rho - gammahat)))^(-1)
    rbet = (1 - mean(X)/(qtp * (k/(n * (1 - tau)))^(gammahat))) * 
      (1/(2 * tau - 1)) * (1 + mopest$beta * gammahat * (1/gammahat - 
                                                           1)^(-mopest$rho) * (1 - tau)^(-mopest$rho)/(gammahat * 
                                                                                                         (1 - mopest$rho - gammahat)))^(-1)
    estimpoint=qtp * (k/(n * (1 - tau)))^gammahat * (1 + ((k/(n * 
                                                                (1 - tau)))^mopest$rho - 1)/mopest$rho * mopest$beta * 
                                                       gammahat * (n/k)^mopest$rho) * (r/rbet)^gammahat * 
      (1 + ((1/gammahat - 1)^(-mopest$rho) * rbet^(-mopest$rho) - 
              1)/mopest$rho * mopest$beta * gammahat * (1 - 
                                                          tau)^(-mopest$rho))/(1 + ((1/gammahat - 1)^(-mopest$rho) * 
                                                                                      r^(-mopest$rho) - 1)/mopest$rho * mopest$beta * gammahat * 
                                                                                 (k/n)^(-mopest$rho))
    psihat=mean((X-qtp)*(X>qtp))
    
    Fbarhat=mean(X>qtp)
               
    psi2hat=min(max(2*Fbarhat*qtp^2*gammahat^2*(1/abs((1-gammahat)*(1-2*gammahat))+Fbarhat^(-mopest$rho)*mopest$beta/mopest$rho*( 1/abs((1-gammahat-mopest$rho)*(1-2*gammahat-mopest$rho)) - 1/abs((1-gammahat)*(1-2*gammahat)) )), psihat^2),sqrt(mean((X-qtp)^4*(X>qtp))))
    
    mhat=mean(X)
    
    alphan=1-Fbarhat
    
    M11phi=k/n*(psi2hat/psihat^2-1)
    
    M12phi=k/n*alphan/(1-alphan)
    
    M22phi=k/n*alphan/(1-alphan)
    
    M11xi=(psihat/qtp)^2*((qtp-mhat)^2*M11phi)/(psihat+Fbarhat*(qtp-mhat))^2
    
    M12xi=gammahat*psihat/qtp*((qtp-mhat)*M12phi)/(psihat+Fbarhat*(qtp-mhat))
    
    M22xi=gammahat^2*M22phi
    
    epsilon=(1/gammahat-1)^(-gammahat)*(1-mhat/qtp)^(-gammahat)*(1-2*k/n)^gammahat*(Fbarhat/(k/n))^gammahat
    
    M11=1/gammahat^2*(epsilon^2*(Fbarhat/(k/n))^2*M11xi-2*epsilon*(Fbarhat/(k/n))*M12xi+M22xi)
    
    M12=1/gammahat*(M12xi-epsilon*(Fbarhat/(k/n))*M11xi)
    
    M22=M11xi
    
    S11=M11*1/(1+Fbarhat/(k/n))^4*(1+8*M11*(1+Fbarhat/(k/n))^(-2)/(k)) #order 2
    
    S12=-1/(1+Fbarhat/(k/n))^2*M12*(1+3*M11*(1+Fbarhat/(k/n))^(-2)/k) #order 2
    
    S22=M22  
    
    nablah1=qtp*(1-2*k/n)*(qtp-mhat)/((mhat-2*qtp*k/n)*1/(1+Fbarhat/(k/n))+2*qtp*k/n-qtp)^2
    
    nablah2=qtp*(1-2*k/n)*(1-1/(1+Fbarhat/(k/n)))*1/(1+Fbarhat/(k/n))*mhat/(2*qtp*k/n*1/(1+Fbarhat/(k/n))+qtp-2*qtp*k/n-1/(1+Fbarhat/(k/n))*mhat)^2
    
    S11prime=nablah1^2*S11+nablah2^2*S22+2*nablah1*nablah2*S12
    
    S12prime=nablah1*S12+nablah2*S22
    
    S22prime=S22
    
    h=gammahat
    
    nabla1=log((2*tau-1)/(1-2*k/n))-log((1-tau)/(k/n))+log(1-mhat/qtp)-log(1-mhat/(qtp*((1-tau)/(k/n))^(-h)))-h*log((1-tau)/(k/n))*mhat/(mhat-qtp*((1-tau)/(k/n))^(-h))
    
    nabla2=1-h*mhat/(qtp*((1-tau)/(k/n))^(-h)-mhat)+h*mhat/(qtp-mhat)
    
    varcorr=nabla1^2*S11prime/log(k/(n*(1-tau)))^2+nabla2^2*S22prime/log(k/(n*(1-tau)))^2+2*nabla1*nabla2*S12prime/log(k/(n*(1-tau)))^2
    
    estimup=estimpoint*exp(-(sqrt(varcorr)/sqrt(k)*log((k/n)/(1-tau))*qnorm((1-ci.level)/2)))
    
    estimdown=estimpoint*exp(-(sqrt(varcorr)/sqrt(k)*log((k/n)/(1-tau))*qnorm((1+ci.level)/2)))
    
    return(list(Lower_bound=estimdown,Point_estimate=estimpoint,Upper_bound=estimup))
    
  }
  if (method == "indirect") {
    gammahat=mopest$EVI[k]
    if(gammahat>1){
      stop("Tail index greater than 1 ! Expectile does probably not exist !")
    }
    qtp = quantile(X, 1 - k/n)*(1/gammahat-1)^(-gammahat)
    r = (1 - mean(X)/qtp) * (n/(n - 2 * k)) * (1 + mopest$beta * 
                                                 gammahat * Fbar(X, qtp)^(-mopest$rho)/(gammahat * (1 - 
                                                                                                      mopest$rho - gammahat)))^(-1)
    rbet = (1 - mean(X)/(qtp * (k/(n * (1 - tau)))^(gammahat))) * 
      (1/(2 * tau - 1)) * (1 + mopest$beta * gammahat * (1/gammahat - 
                                                           1)^(-mopest$rho) * (1 - tau)^(-mopest$rho)/(gammahat * 
                                                                                                         (1 - mopest$rho - gammahat)))^(-1)
    estimpoint=(1/gammahat - 1)^(-gammahat) * quantile(X, 1 - 
                                                         k/n) * (k/(n * (1 - tau)))^gammahat * (1 + ((k/(n * 
                                                                                                           (1 - tau)))^mopest$rho - 1)/mopest$rho * mopest$beta * 
                                                                                                  gammahat * (n/k)^mopest$rho) * (1 + ((1/gammahat - 
                                                                                                                                          1)^(-mopest$rho) * rbet^(-mopest$rho) - 1)/mopest$rho * 
                                                                                                                                    mopest$beta * gammahat * (1 - tau)^(-mopest$rho))/(rbet^gammahat)
    quanthat=quantile(X,1-k/n)
    
    expectint=expect(X,1-k/n)
    
    mhat=mean(X)
    
    V11=1
    
    V12=(1/(1-gammahat)-log(1/gammahat-1)) #ordre 1
    
    V12=V12+(3*gammahat-1)/(2*(1-gammahat)^3*k) #order 2
    
    V12=V12+3*(10*gammahat^3-10*gammahat^2+5*gammahat-1)/(4*(1-gammahat)^5*k^2) #order 3
    
    V22=(1+(1/(1-gammahat)-log(1/gammahat-1))^2) #order 1
    
    V22=V22+(1/(1-gammahat)-log(1/gammahat-1))*(3*gammahat-1)/((1-gammahat)^3*k)+0.5/((1-gammahat)^4*k) #order 2
    
    V22=V22+36*(10*gammahat^3-10*gammahat^2+5*gammahat-1)*(1-(1-gammahat)*log(1/gammahat-1))/(24*(1-gammahat)^6*k^2)+(6*gammahat^2-4*gammahat+1)/((1-gammahat)^6*k^2)+20*(3*gammahat-1)^2/(48*(1-gammahat)^6*k^2) #order 3
    
    nabla1=1-gammahat*mhat/((n*(1-tau)/k)^(-gammahat)*(1/gammahat-1)^(-gammahat)*quanthat-mhat)+(log(2*tau-1)-log(1-mhat/((n*(1-tau)/k)^(-gammahat)*(1/gammahat-1)^(-gammahat)*quanthat)))/log(k/(n*(1-tau)))
    
    nabla2=(1-mhat*gammahat/(quanthat*(1/gammahat-1)^(-gammahat)*(k/(n*(1-tau)))^gammahat-mhat))/log(k/(n*(1-tau)))
    
    varcorr=nabla1^2*V11*gammahat^2+nabla2^2*V22*gammahat^2+2*nabla1*nabla2*V12*gammahat^2
    
    estimup=estimpoint*exp(-(sqrt(varcorr)/sqrt(k)*log((k/n)/(1-tau))*qnorm((1-ci.level)/2)))
    
    estimdown=estimpoint*exp(-(sqrt(varcorr)/sqrt(k)*log((k/n)/(1-tau))*qnorm((1+ci.level)/2)))
    
    return(list(Lower_bound=estimdown,Point_estimate=estimpoint,Upper_bound=estimup))
  }
}
