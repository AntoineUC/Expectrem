logspacqqplot=function(X,k=trunc(length(X)/10),grid=seq(0.01,0.99,0.01),weighted=T,add.line=F){
  
  n=length(X)
  
  if(!is.logical(weighted)){
    stop("weighted must be boolean.")
  }
  if(!is.logical(add.line)){
    stop("add.line must be boolean.")
  }
  if (length(k) > 1) {
    stop("k must be of length 1.")
  }
  if (k > n - 1 || k < 1) {
    stop("k must be between 1 and n-1.")
  }
  
  if(any(grid>=1) || any(grid<=0)){
    stop("grid must be strictly between 0 and 1")
  }
  
  if(weighted==T){
    plot(y=quantile((1:k)*log(quantile(X,1-(1:k-1)/n,type=1)/quantile(X,1-1:k/n,type=1)),grid),x=(-log(1-grid)),pch=16,main="QQ-plot",ylab="Log-spacings quantiles",xlab="Exponential quantiles")
  }
  if(weighted==F){
    plot(y=quantile(log(quantile(X,1-1:k/n,type=1)/quantile(X,1-k/n,type=1)),grid),x=(-log(1-grid)),pch=16,main="QQ-plot",ylab="Log-spacings quantiles",xlab="Exponential quantiles")
  }
  if(add.line==T){
    abline(a=0,b=mop(X,k=k,p=0,method="RBMOP")$EVI)
  }
}
