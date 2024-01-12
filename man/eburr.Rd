\name{eburr}
\alias{eburr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of the Burr distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the Burr distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eburr(probs, a, b, start.pt=((1-probs)^(-1/b)-1)^(1/a), tol=1e-08, maxiter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{a}{
Parameter of the Burr distribution. Must be strictly positive.
}
\item{b}{
Parameter of the Burr distribution. Must be strictly greater than 1/a.
}
\item{start.pt}{
Vector of starting points. Default: quantiles of the Burr distribution.
}
\item{tol}{
Tolerance in the stopping criterion. Default: 1e-08.
}
\item{maxiter}{
Maximum number of iterations. Default: 100.
}
}
\details{
The considered distribution has c.s.f. (1+ x^a)^(-b). Note that this distribution is sometimes parametrized with gamma=1/(a*b) and rho=-1/b.%%  ~~ If necessary, more details than the description above ~~
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
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~Burr }% __ONLY ONE__ keyword per line
