logspaceeplot=function(X,k=trunc(length(X)/10),grid=seq(0.01,0.99,0.01),add.line=F){
  
  n=length(X)
  
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
  
  plot(y=expect(log(quantile(X,1-1:k/n,type=1)/quantile(X,1-k/n,type=1)),grid),x=eexp(grid),pch=16,main="EE-plot",ylab="Log-spacings expectiles",xlab="Exponential expectiles")
  if(add.line==T){
    abline(a=0,b=mop(X,k=k,p=0,method="RBMOP")$EVI)
  }
}