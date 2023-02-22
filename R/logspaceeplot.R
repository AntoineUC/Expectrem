logspaceeplot=function (X, k = trunc(length(X)/10), weighted = T, add.line = F) 
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
    plot(y = expect((1:k) * log(quantile(X, 1 - (1:k - 1)/n, 
                                         type = 1)/quantile(X, 1 - 1:k/n, type = 1)), (1:k)/(k+1)), 
         x = eexp((1:k)/(k+1)), pch = 16, main = "EE-plot", ylab = "Log-spacings expectiles", 
         xlab = "Exponential expectiles")
  }
  if (weighted == F) {
    plot(y = expect(log(quantile(X, 1 - 1:k/n, type = 1)/quantile(X, 
                                                                  1 - k/n, type = 1)), (1:k)/(k+1)), x = eexp((1:k)/(k+1)), pch = 16, 
         main = "EE-plot", ylab = "Log-spacings expectiles", 
         xlab = "Exponential expectiles")
  }
  if (add.line == T) {
    abline(a = 0, b = mop(X, k = k, p = 0, method = "RBMOP")$EVI)
  }
}
