\name{extES}
\alias{extES}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extreme Expected Shortfall estimation
%%  ~~function to do ... ~~
}
\description{This function computes several estimators for extreme Expected Shortfalls.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extES(X,k,tau,estim="Moment",method="direct",ci.level=0.95,n.bootstrap=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
A vector of number of upper order statistics.
}
  \item{tau}{
Numeric (extreme) probability with value in [0,1].
}
  \item{estim}{
ES estimator used : Moment estimator ("Moment", default) tailored for all domains of attraction, or Weissman/Hill estimator ("Hill", default), only suitable for heavy-tailed distribution (positive tail index).
}
  \item{method}{
Direct ("direct", default) or indirect ("indirect") estimator.
}
  \item{ci.level}{
Confidence level. Default: 0.95.
}
  \item{n.bootstrap}{
Size of the bootstrap sample simulated (only if method="direct" and estim="Moment", default 10,000).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Lower_bound }{Lower bound of the confidence interval.}
  \item{Point_estimate }{Expected Shortfall estimate.}
   \item{Upper_bound }{Upper bound of the confidence interval.}
}
\references{
Daouia, A., Stupfler, G. & Usseglio-Carleve, A. (2024). A unified framework for extreme Expected Shortfall inference.
}
\author{
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Expected Shortfall }% use one of  RShowDoc("KEYWORDS")
