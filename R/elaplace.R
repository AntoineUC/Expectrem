elaplace=function(probs,location=0,scale=1){
  
  if(any(probs>=1) || any(probs<=0)){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if(scale<=0){
    stop("scale must be strictly positive.")
  }
  
  return(location+scale*sign(2*probs-1)*lambertWp(abs(2*probs-1)/(1-abs(2*probs-1))))
}