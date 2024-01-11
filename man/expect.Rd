\name{expect}
\alias{expect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Empirical expectiles
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes empirical expectiles corresponding to the given probabilities.
}
\usage{
expect(X, probs, tol = 1e-08, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Data vector.
}
  \item{probs}{
%%     ~~Describe \code{alphas} here~~
Numeric vector of probabilities with values in (0,1).
}
\item{tol}{
Tolerance (default: 1e-08).
}
\item{maxiter}{
Maximum number of iterations (default: 100).
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
\code{\link{quantile}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
