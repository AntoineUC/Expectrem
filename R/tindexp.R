tindexp=function(X,k){
  n=length(X)
  return(1/(1+Fbar(X,expect(X,1-k/n))/(k/n)))
}
