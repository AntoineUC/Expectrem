############################################################
####                Expectile functions                 ####
############################################################


library(Rcpp)

library(RcppArmadillo)

library(evt0)

#The following function computes the expectile empirical c.s.f.

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            double F2(arma::vec& X, const double y) {
            return(sum((X-y)%(X>y))/sum(abs(X-y)));
            }
            ") #expectile c.s.f.

#The following function computes the empirical c.s.f.

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec Fbar(arma::vec& X, arma::vec& y) {
            arma::mat M=repmat(X,1,y.n_elem);
            arma::mat Y=repmat(trans(y),X.n_elem,1);
            arma::vec N(y.n_elem);
            N.ones();
            double n=X.n_elem;
            N=n*N;
            return(sum(trans(M>Y),1)/N);
            }
            ") #quantile c.s.f.

#The following function computes the empirical expectiles

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec expect(arma::vec& X, arma::vec& alphas) {
            arma::mat M=repmat(X,1,alphas.n_elem);
            arma::vec a(arma::size(alphas));
            arma::vec b(arma::size(alphas));
            a.zeros();
            b.zeros();
            arma::vec gap(arma::size(alphas));
            gap.ones();
            arma::mat A=repmat(trans(a),X.n_elem,1);
            while(gap.max()>0.00001){
              A=repmat(trans(a),X.n_elem,1);
              b=mean(trans(M),1)-alphas%mean(trans(M-A)%trans(M<A),1)-(1-alphas)%mean(trans(M-A)%trans(M>A),1);
              gap=abs(b-a);
              a=b;
            }
            return(b);
            }
            ") #empirical expectile function

#The following function computes the expectile based tail index estimator

tindexp=function(X,k){
  n=length(X)
  return(1/(1+Fbar(X,expect(X,1-k/n))/(k/n)))
}

#The following function computes the bias reduced version of the expectile based tail index estimator

tindexpbr=function(X,k){
  
  n=length(X)
  
  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")
  
  gammhill=mopest$EVI[k]
  
  qtp=expect(X,1-k/n)
  
  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammhill*Fbar(X,qtp)^(-mopest$rho)/(gammhill*(1-mopest$rho-gammhill)))^(-1)
  
  return(1/(1+Fbar(X,qtp)*n/k*1/r))
}

#The following function computes the Weissman estimator for extreme quantiles using different tail index estimators

extQuant=function(X,k,tau,estim="Hill"){
  n=length(X)
  if(estim=="Hillbr"){
  gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")$EVI[k]
  }
  if(estim=="Hill"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k)
  }
  if(estim=="tindexpbr"){
    gammahat=tindexpbr(X,k)
  }
  return(quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat)
}

#The following function computes a bias reduced version of the Weissman estimator for extreme quantiles using different tail index estimators

extQuantbr=function(X,k,tau,estim="Hill"){
  n=length(X)
  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")
  if(estim=="Hillbr"){
    gammahat=mopest$EVI[k]
  }
  if(estim=="Hill"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k)
  }
  if(estim=="tindexpbr"){
    gammahat=tindexpbr(X,k)
  }
  return(quantile(X,1-k/n)*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho))
}

#The following function computes the Weissman type estimator for extreme expectiles using different tail index estimators

extExpect=function(X,k,tau,estim="Hill"){
  n=length(X)
  if(estim=="Hillbr"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")$EVI[k]
  }
  if(estim=="Hill"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k)
  }
  if(estim=="tindexpbr"){
    gammahat=tindexpbr(X,k)
  }
  return(expect(X,1-k/n)*(k/(n*(1-tau)))^gammahat)
}

#The following function computes a bias reduced version of the Weissman type estimator for extreme expectiles using different tail index estimators

extExpectbr=function(X,k,tau,estim="Hill"){
  n=length(X)
  mopest=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="RBMOP")
  if(estim=="Hillbr"){
    gammahat=mopest$EVI[k]
  }
  if(estim=="Hill"){
    gammahat=mop(X[which(X>0)], 1:(length(which(X>0))-1), 0, method ="MOP")$EVI[k]
  }
  if(estim=="tindexp"){
    gammahat=tindexp(X,k)
  }
  if(estim=="tindexpbr"){
    gammahat=tindexpbr(X,k)
  }
  
  qtp=expect(X,1-k/n)
  
  r=(1-mean(X)/qtp)*(n/(n-2*k))*(1+mopest$beta*gammahat*Fbar(X,qtp)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)
  
  rbet=(1-mean(X)/(qtp*(k/(n*(1-tau)))^(gammahat)))*(1/(2*tau-1))*(1+mopest$beta*gammahat*(1/gammahat-1)^(-mopest$rho)*(1-tau)^(-mopest$rho)/(gammahat*(1-mopest$rho-gammahat)))^(-1)
  
  return(qtp*(k/(n*(1-tau)))^gammahat*(1+((k/(n*(1-tau)))^mopest$rho-1)/mopest$rho*mopest$beta*gammahat*(n/k)^mopest$rho)*(r/rbet)^gammahat*(1+((1/gammahat-1)^(-mopest$rho)*rbet^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(1-tau)^(-mopest$rho))/(1+((1/gammahat-1)^(-mopest$rho)*r^(-mopest$rho)-1)/mopest$rho*mopest$beta*gammahat*(k/n)^(-mopest$rho)))
}
