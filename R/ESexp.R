ESexp=function(probs,lambda=1){
  
  if(any(probs>=1) || any(probs<=0)){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(length(lambda)>1){
    stop("lambda must be of size 1.")
  }
  
  if(lambda<=0){
    stop("lambda must be strictly positive.")
  }
  
  return((1-log(1-probs))/lambda)
}