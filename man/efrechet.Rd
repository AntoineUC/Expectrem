\name{efrechet}
\alias{efrechet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of the Frechet distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the expectiles of the Frechet distribution.
}
\usage{
efrechet(probs,a=2,niter=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
%%     ~~Describe \code{x} here~~
Numeric vector of probabilities with values in (0,1).
}
\item{a}{
%%     ~~Describe \code{x} here~~
Parameter (strictly greater than 1) of the Frechet distribution.
}
\item{niter}{
Number of Newton-Raphson iterations. Default: 50.
}
}
\details{
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
%%  ~~who you are~~
Antoine Usseglio-Carleve
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.