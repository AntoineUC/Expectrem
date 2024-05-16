\name{momentindex}
\alias{momentindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Moment estimator of the tail index
%%  ~~function to do ... ~~
}
\description{This function computes the Moment tail index estimator introduced by Dekkers.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
momentindex(x,k,ci.level=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data vector.
}
\item{k}{
A vector of number of upper order statistics.
}
\item{ci.level}{
Confidence level. Default: 0.95.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Lower_bound }{Lower bound of the confidence interval.}
  \item{Point_estimate }{Tail index estimate.}
   \item{Upper_bound }{Upper bound of the confidence interval.}
}
\references{
Dekkers, A. L., Einmahl, J. H., & De Haan, L. (1989). A moment estimator for the index of an extreme-value distribution. The Annals of Statistics, 1833-1855.
}
\author{Antoine Usseglio-Carleve
%%  ~~who you are~~
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
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
