\name{epareto}
\alias{epareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of Pareto distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the Pareto distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
epareto(probs, alpha, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{alpha}{
Shape parameter strictly greater than 1.
}
\item{start.pt}{
  Choice of the starting points. Mean of the Pareto distribution if "mean", quantiles (levels probs) of the Pareto distribution if "quantile". The starting point is chosen manually if "custom".  
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
The considered distribution has c.s.f. y^(-alpha) for y>1.%%  ~~ If necessary, more details than the description above ~~
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

\keyword{ ~expectiles }
\keyword{ ~pareto }
