\name{extExpect}
\alias{extExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extreme expectile estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the Weissman-type estimator for extreme expectiles.
}
\usage{
extExpect(X, k="kopt", tau, estim = "Hill", method="direct", br = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Data vector.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
A vector of number of upper order statistics. Default: optimal k given in Girard et al. (2020).
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
Numeric (extreme) probability with value in [0,1].
}
  \item{estim}{
%%     ~~Describe \code{estim} here~~
Tail index estimator used : Hill estimator ("Hill", default) or proportionality estimator ("tindexp").
}
  \item{method}{
%%     ~~Describe \code{estim} here~~
Direct or indirect (quantile based) estimator.
}
\item{br}{
Bias reduced version if TRUE. Default: FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Daouia, A., Girard, S., and Stupfler, G.  (2018). Estimation of tail risk based on extreme expectiles. Journal of the Royal Statistical Society: Series B, 80(2): 263-292.
Girard, S., Stupfler, G., and Usseglio-Carleve, A. (2020). On second-order automatic bias reduction for extreme expectile estimation, preprint.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~extreme }% __ONLY ONE__ keyword per line
