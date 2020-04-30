\name{extExpectbr}
\alias{extExpectbr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bias-reduced extreme expectile estimation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes a bias-reduced version of the Weissman-type estimator for extreme expectiles.
}
\usage{
extExpectbr(X, k, tau, estim = "Hill")
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
Daouia, A., Girard, S., and Stupfler, G.  (2018). Estimation of tail risk based on extreme expectiles. Journal of the Royal Statistical Society: Series B, 80(2): 263-292.
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
\code{\link{extExpect}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~extreme }% __ONLY ONE__ keyword per line
\keyword{ ~reduced-bias }% __ONLY ONE__ keyword per line
