\name{tindexp}
\alias{tindexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Proportionality tail index estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes the proportionality estimator of the extreme value index.
}
\usage{
tindexp(X, k, br = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
Data vector.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
A vector of number of upper order statistics.
}
\item{br}{
%%     ~~Describe \code{k} here~~
Bias-reduced version if TRUE. Default: FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This estimator is defined as 1/(1+Fbar(X,e(X,1-k/n))/(k/n)), where Fbar is the empirical c.s.f., and e(X,1-k/n) the empirical expectile of level 1-k/n.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A vector of k estimates of the tail index.
}
\references{
%% ~put references to the literature/web site here ~
Girard, S., Stupfler, G., and Usseglio-Carleve, A. (2019). Nonparametric extreme conditional expectile estimation, preprint. https://hal.archives-ouvertes.fr/hal-02114255.
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
\keyword{ ~hill }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~expectile }% __ONLY ONE__ keyword per line
