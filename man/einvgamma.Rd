\name{einvgamma}
\alias{einvgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of Inverse-Gamma distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the Inverse-Gamma distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
einvgamma(probs, alpha, lambda, tol=1e-08, maxiter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{alpha}{
Shape parameter strictly greater than 1. 
}
\item{lambda}{
Scale parameter strictly greater than 0. 
}
 \item{tol}{
Tolerance in the stopping criterion. Default: 1e-08.
}
\item{maxiter}{
Maximum number of iterations. Default: 100.
}
}
\details{
The considered distribution has p.d.f. lambda*alpha*y^(-alpha-1)*exp(-lambda/y)/Gamma(alpha) for y>0.%%  ~~ If necessary, more details than the description above ~~
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
