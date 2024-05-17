\name{ESweibull}
\alias{ESweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected Shortfalls of Weibull distribution
  %%  ~~function to do ... ~~
}
\description{
  This function computes the Expected Shortfalls of the Weibull distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  ESweibull(probs, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
    Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
  }
  \item{k}{
    Shape parameter strictly positive. 
  }
}
\details{
  The considered distribution has c.s.f. exp(-y^k) for y>0.%%  ~~ If necessary, more details than the description above ~~
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
\keyword{ ~weibull }