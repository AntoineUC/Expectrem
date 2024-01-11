\name{echisq}
\alias{echisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of chi-squared distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the chi-squared distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
echisq(probs, df, tol=1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{df}{
Degrees of freedom. Must be strictly positive.
}
 \item{tol}{
Tolerance in the stopping criterion. Default: 1e-08.
}
\item{maxiter}{
Maximum number of iterations. Default: 100.
}
}
\details{
The considered distribution has p.d.f. y^(df/2-1)*exp(-y/2)/(2^(df/2)*Gamma(df/2)) for y>0.%%  ~~ If necessary, more details than the description above ~~
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
