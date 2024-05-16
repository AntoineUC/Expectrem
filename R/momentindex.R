momentindex=function(X,k,ci.level=0.95){
  
  if (length(ci.level) > 1) {
    stop("ci.level must be of length 1.")
  }
  
  if (ci.level >= 1 || ci.level <= 0) {
    stop("ci.level must be strictly between 0 and 1.")
  }
  
  n=length(X)
  
  if(any(k>n-1) || any(k<1)){
    stop("k must be between 1 and n-1.")
  }
  
  fM1=function(k){
    return(mean(log(quantile(X,1-(0:(k-1))/n,type=1)/quantile(X,1-k/n,type=1))))
  }
  
  fM2=function(k){
    return(mean(log(quantile(X,1-(0:(k-1))/n,type=1)/quantile(X,1-k/n,type=1))^2))
  }
  
  M1=sapply(k,fM1)
  
  M2=sapply(k,fM2)
  
  M=M1+1-0.5*(1-M1^2/M2)^(-1)
  
  variance=(M^2+1)*(M>=0)+((1-M)^2*(1-2*M)*(1-M+6*M^2))/((1-3*M*(M<0))*(1-4*M*(M<0)))*(M<0)
  
  Mdown=M-qnorm(1-(1-ci.level)/2)*sqrt(variance/k)
    
  Mup=M+qnorm(1-(1-ci.level)/2)*sqrt(variance/k)
  
  return(list(Lower_bound=Mdown,Point_estimate=M,Upper_bound=Mup))  
  
}