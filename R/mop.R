mop=function (x, k, p, method = c("MOP", "RBMOP")) 
{
  n = length(x)
  if (n < 2) {
    stop("Data vector x with at least two sample points required.")
  }
  if (is.null(p) || any(is.na(p))) {
    stop("p is not specified")
  }
  if (is.null(k) || any(is.na(k))) {
    stop("k is not specified")
  }
  if (any(k < 1) || any(k > n) || any(k == n) || !is.numeric(k) || 
      k != as.integer(k)) {
    stop("Each k must be integer and greater than or equal to 1 and less than sample size.")
  }
  if (!is.numeric(x)) {
    stop("Some of the data points are not real.")
  }
  if (!is.numeric(p)) {
    stop("p must be a real.")
  }
  osx <- sort(x[x > 0])
  method = match.arg(method)
  if (method == "MOP") {
    EVI.est = mop.MOP(osx, k, p)
  }
  else if (method == "RBMOP") {
    EVI.est = mop.RBMOP(osx, k, p)
  }
  colnames(EVI.est$EVI) = p
  rownames(EVI.est$EVI) = k
  return(EVI.est)
}