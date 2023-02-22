\name{logspacqqplot}
\alias{logspacqqplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exponential Quantile-Quantile plot of the log-spacings
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logspacqqplot(X,k=trunc(length(X)/10),weighted=T,add.line=F)
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
  \item{weighted}{
%%     ~~Describe \code{x} here~~
Boolean. If TRUE (default), the weighted log-spacings i*log(X_{n-i+1,n}/X_{n-i,n}) for i between 1 and k are computed. If FALSE, the log-spacings log(X_{n-i,n}/X_{n-k,n}) are considered.
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
\code{\link{logspaceeplot}}
}
\examples{

}

