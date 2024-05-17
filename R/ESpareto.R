ESpareto=function (probs, alpha) 
{
  
  if (min(probs) <= 0 || max(probs) >= 1) {
    stop("only asymmetries between 0 and 1 allowed.")
  }
  
  if (length(alpha) > 1) {
    stop("alpha must be of length 1.")
  }
  
  if (alpha <= 1) {
    stop("alpha must be greater than 1.")
  }

  return(alpha/(alpha-1)*(1-probs)^(-1/alpha))
}