\name{et}
\alias{et}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of Student distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the theoretical expectiles of the Student distribution.
}
\usage{
et(probs, df = 2, mu = 0, sigma = 1, niter = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
%%     ~~Describe \code{probs} here~~
Numeric vector of probabilities with values in [0,1].
}
\item{df}{
%%     ~~Describe \code{probs} here~~
Degrees of freedom parameter. Default: 2.
}
\item{mu}{
%%     ~~Describe \code{niter} here~~
Location parameter. Default: 0.
}
\item{sigma}{
%%     ~~Describe \code{niter} here~~
Scale parameter. Default: 1.
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
Number of iterations of the fixed point algorithm. Default: 20.
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
Maume-Deschamps, V., Rulliere, D., and Usseglio-Carleve, A. (2018). Spatial expectile predictions for elliptical random fields. Methodology and Computing in Applied Probability, 20(2):643-671.
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
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~student }% __ONLY ONE__ keyword per line
