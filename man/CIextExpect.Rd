\name{CIextExpect}
\alias{CIextExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extreme expectile estimator with confidence interval
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the (direct and indirect) Weissman type estimators of extreme expectiles, and returns the associated confidence intervals.
}
\usage{
CIextExpect(X, k = trunc(length(X)/10), tau, method = "direct",ci.level=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
Data vector.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
A vector of number of upper order statistics. Default: 10\% of the data.
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
Numeric (extreme) probability with value in [0,1].
}
  \item{method}{
%%     ~~Describe \code{estim} here~~
Direct or indirect (quantile based) estimator. Default: direct.
}
  \item{ci.level}{
%%     ~~Describe \code{estim} here~~
Confidence level. Default: 0.95.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Lower_bound }{Lower bound of the confidence interval.}
  \item{Point_estimate }{Extreme expectile estimate.}
   \item{Upper_bound }{Upper bound of the confidence interval.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Daouia, A., Stupfler, G. and Usseglio-Carleve, A. (2023). Bias- and variance-corrected asymptotic Gaussian inference about extreme expectiles.
}
\author{
%%  ~~who you are~~
Antoine Usseglio-Carleve
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
