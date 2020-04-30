\name{extQuantbr}
\alias{extQuantbr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bias-reduced extreme quantile estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes a bias-reduced version of the Weissman-type estimator for extreme quantiles.
}
\usage{
extQuantbr(X, k, tau, estim = "Hill")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Data vector.
}
  \item{k}{
%%     ~~Describe \code{k} here
A vector of number of upper order statistics.
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
Numeric (extreme) probability with value in [0,1].
}
  \item{estim}{
%%     ~~Describe \code{estim} here~~
Tail index estimator used : classical Hill estimator ("Hill", default), bias-reduced Hill estimator ("Hillbr"), proportionality estimator ("tindexp") and bias-reduced proportionality estimator ("tindexpbr").
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
Gomes, M. I. and Pestana, D. (2007). A sturdy reduced-bias extreme quantile (VaR) estimator. Journal of the American Statistical Association, 102(477): 280-292.

Weissman, I. (1978). Estimation of parameters and large quantiles based on the k largest observations. Journal of the American Statistical Association, 73(364): 812-815.
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
\code{\link{extQuant}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~quantile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~extreme }% __ONLY ONE__ keyword per line
\keyword{ ~reduced-bias }% __ONLY ONE__ keyword per line
