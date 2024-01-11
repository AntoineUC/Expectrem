\name{egpd}
\alias{egpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of Generalized Pareto Distribution
%%  ~~function to do ... ~~
}
\description{
This function computes expectiles of the Generalized Pareto Distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
egpd(probs, alpha, tol=1e-08, maxiter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{alpha}{
Parameter strictly greater than 1.
}
 \item{tol}{
Tolerance in the stopping criterion. Default: 1e-08.
}
\item{maxiter}{
Maximum number of iterations. Default: 100.
}
}
\details{
The considered GPD has c.s.f. (1+ y/alpha)^(-alpha). Note that the GPD is heavy-tailed with first and second order tail indices gamma=1/alpha and rho=-1/alpha, respectively.%%  ~~ If necessary, more details than the description above ~~
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
\keyword{ ~expectiles }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~gpd }% __ONLY ONE__ keyword per line
