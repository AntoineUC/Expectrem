\name{ESkumaraswamy}
\alias{ESkumaraswamy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected Shortfalls of the Kumaraswamy distribution
  %%  ~~function to do ... ~~
}
\description{
  This function computes the Expected Shortfalls of the Kumaraswamy distribution.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  ESkumaraswamy(probs, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
    Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
  }
  \item{a}{
    First parameter. Must be strictly positive.
  }
  \item{b}{
    Second parameter. Must be strictly positive.
  }
}
\details{
  The considered distribution has c.d.f. 1-(1-y^a)^b for y in [0,1].
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

\keyword{ ~Expected Shortfall }