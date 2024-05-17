ESfrechet=function(probs, a){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(a) > 1) {
    stop("a must be of length 1.")
  }
  
  if(a<=1){
    stop("Expectiles do not exist if a is less than 1.")
  }
  
  return(pgamma(-log(probs),1-1/a)*gamma(1-1/a)/(1-probs))
}