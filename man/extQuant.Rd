\name{extQuant}
\alias{extQuant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extreme quantile estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the Weissman-type estimator for extreme quantiles.
}
\usage{
extQuant(X, k, tau, estim = "Hill", br = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Data vector.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
A vector of number of upper order statistics.
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
Numeric (extreme) probability with value in [0,1].
}
  \item{estim}{
%%     ~~Describe \code{estim} here~~
Tail index estimator used : Hill estimator ("Hill", default) or proportionality estimator ("tindexp").
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
Weissman, I. (1978). Estimation of parameters and large quantiles based on the k largest observations. Journal of the American Statistical Association, 73(364):812-815.
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
\keyword{ ~quantile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~extreme }% __ONLY ONE__ keyword per line
