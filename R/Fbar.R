Fbar=function(X,y){
  return(apply((t(matrix(rep(X,length(y)),length(X),length(y)))-matrix(rep(y,length(X)),length(y),length(X))>0),1,mean))
}
