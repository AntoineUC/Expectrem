\name{ebeta}
\alias{ebeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of Beta distribution
  %%  ~~function to do ... ~~
}
\description{
  This function computes the expectiles of the Beta distribution.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  ebeta(probs, shape1, shape2, start.pt="quantile", tol=1e-08, maxiter=100, x0=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
    Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
  }
  \item{shape1}{
    First shape parameter. Must be strictly positive.
  }
  \item{shape2}{
    Second shape parameter. Must be strictly positive.
  }
  \item{start.pt}{
    Choice of the starting points. Mean of the Beta distribution if "mean", quantiles (levels probs) of the Beta distribution if "quantile". The starting point is chosen manually if "custom".  
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
  The considered distribution has p.d.f. y^(shape1-1)*(1-y)^(shape2-1)/B(shape1,shape2) for y in [0,1].
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

\keyword{ ~expectiles }
