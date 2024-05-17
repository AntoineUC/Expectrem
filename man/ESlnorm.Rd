\name{ESlnorm}
\alias{ESlnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected Shortfalls of Lognormal distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the theoretical Expected Shortfalls of the Lognormal distribution.
}
\usage{
ESlnorm(probs, mu = 0, sigma = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
%%     ~~Describe \code{probs} here~~
Numeric vector of probabilities with values in [0,1].
}
\item{mu}{
%%     ~~Describe \code{niter} here~~
Meanlog parameter. Default: 0.
}
\item{sigma}{
%%     ~~Describe \code{niter} here~~
sdlog parameter. Default: 1.
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
Antoine Usseglio-Carleve.
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
\keyword{ ~Expected Shortfall }% use one of  RShowDoc("KEYWORDS")