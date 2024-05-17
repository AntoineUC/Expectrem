ESt=function(probs, df, mu = 0, sigma = 1){
  
  if (min(probs) <= 0 || max(probs) >= 1){
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(df) > 1) {
    stop("df must be of length 1.")
  }
  
  if (length(mu) > 1) {
    stop("mu must be of length 1.")
  }
  
  if (length(sigma) > 1) {
    stop("sigma must be of length 1.")
  }
  
  if(sigma<0){
    stop("sigma must be positive.")
  }
  
  if(df<=1){
    stop("df must greater than 1.")
  }
  
  return(mu+sigma*gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2))*df/(df-1)*(1+qt(probs,df)^2/df)^((1-df)/2)*1/(1-probs))
}