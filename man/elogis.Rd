\name{elogis}
\alias{elogis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of Logistic distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the theoretical expectiles of the Logistic distribution.
}
\usage{
elogis(probs, mu = 0, sigma = 1, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
%%     ~~Describe \code{x} here~~
Numeric vector of probabilities with values in [0,1].
}
\item{mu}{
%%     ~~Describe \code{x} here~~
Location parameter. Default: 0.
}
\item{sigma}{
%%     ~~Describe \code{x} here~~
Scale parameter. Default: 1.
}
\item{start.pt}{
  Choice of the starting points. 0 if "mean", quantiles (levels probs) of the Logistic distribution if "quantile". The starting point is chosen manually if "custom".  
}
\item{tol}{
%%     ~~Describe \code{x} here~~
Tolerance in the stopping criterion. Default: 1e-08.
}
\item{maxiter}{
Maximum number of iterations. Default: 100.
}
\item{x0}{
  Value of the starting point (only if start.pt="custom").
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
Maume-Deschamps, V., Rulliere, D., and Usseglio-Carleve, A. (2018). Spatial expectile predictions for elliptical random fields. Methodology and Computing in Applied Probability, 20(2):643-671.
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
\keyword{ ~logistic }% __ONLY ONE__ keyword per line
