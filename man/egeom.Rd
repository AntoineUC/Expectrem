\name{egeom}
\alias{egeom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of geometric distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the geometric distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
egeom(probs,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{p}{
Success probability. Must be between 0 and 1.
}
}
\details{
The geometric distribution with prob = p has density p(x) = p (1-p)^x for x = 0, 1, 2, …, 0 < p ≤ 1.%%  ~~ If necessary, more details than the description above ~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~geometric }% __ONLY ONE__ keyword per line
