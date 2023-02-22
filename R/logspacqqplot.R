logspacqqplot=function (X, k = trunc(length(X)/10), weighted = T, add.line = F) 
{
  n = length(X)
  if (!is.logical(weighted)) {
    stop("weighted must be boolean.")
  }
  if (!is.logical(add.line)) {
    stop("add.line must be boolean.")
  }
  if (length(k) > 1) {
    stop("k must be of length 1.")
  }
  if (k > n - 1 || k < 1) {
    stop("k must be between 1 and n-1.")
  }
  if (weighted == T) {
    plot(y = quantile((1:k) * log(quantile(X, 1 - (1:k - 
                                                     1)/n, type = 1)/quantile(X, 1 - (1:k)/n, type = 1)), 
                      (1:k)/(k+1)), x = (-log(1 - (1:k)/(k+1))), pch = 16, main = "QQ-plot", 
         ylab = "Log-spacings quantiles", xlab = "Exponential quantiles")
  }
  if (weighted == F) {
    plot(y = quantile(log(quantile(X, 1 - 1:k/n, type = 1)/quantile(X, 
                                                                    1 - k/n, type = 1)), (1:k)/(k+1)), x = (-log(1 - (1:k)/(k+1))), 
         pch = 16, main = "QQ-plot", ylab = "Log-spacings quantiles", 
         xlab = "Exponential quantiles")
  }
  if (add.line == T) {
    abline(a = 0, b = mop(X, k = k, p = 0, method = "RBMOP")$EVI)
  }
}
