\name{edagum}
\alias{edagum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of the Dagum distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the Dagum distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
edagum(probs, a, b, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{a}{
Parameter of the Dagum distribution. Must be strictly greater than 1.
}
\item{b}{
Parameter of the Dagum distribution. Must be strictly positive.
}
\item{start.pt}{
  Choice of the starting points. Mean of the Dagum distribution if "mean", quantiles (levels probs) of the Dagum distribution if "quantile". The starting point is chosen manually if "custom".  
}
\item{tol}{
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
The considered distribution has c.d.f. (1+ y^(-a))^(-b). Note that this distribution is heavy-tailed with first and second tail indices gamma=1/a and rho=-1, respectively.%%  ~~ If necessary, more details than the description above ~~
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
}
\author{
Antoine Usseglio-Carleve%%  ~~who you are~~
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
