\name{ESpareto}
\alias{ESpareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected Shortfalls of Pareto distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the Expected Shortfalls of the Pareto distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ESpareto(probs, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{alpha}{
Shape parameter strictly greater than 1.
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

\keyword{ ~Expected Shortfall }
\keyword{ ~pareto }