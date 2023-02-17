\name{logspaceeplot}
\alias{logspaceeplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponential Expectile-Expectile plot of the log-spacings
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function produces an exponential EE plot of the values of the log spacings of the k greatest values of X.
}
\usage{
logspaceeplot(X,k=trunc(length(X)/10),grid=seq(0.01,0.99,0.01),add.line=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
Data vector.
}
  \item{k}{
%%     ~~Describe \code{x} here~~
Number of upper order statistics. Default: 10\% of the data.
}
  \item{grid}{
%%     ~~Describe \code{x} here~~
A grid over (0,1) where expectiles are estimated. Default: seq(0.01,0.99,0.01).
}
  \item{add.line}{
%%     ~~Describe \code{x} here~~
Add a straight line with slope the estimated tail index (bias-reduced Hill estimator) if TRUE. Default: FALSE.
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
