\name{eburr}
\alias{eburr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of the Burr distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the theoretical expectiles of the Burr distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eburr(probs, gamma=0.5,rho=-1,method="mc")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{gamma}{
First order parameter in (0,1). Default: 0.5.
}
\item{rho}{
Second order parameter in (-Inf,0]. Default: -1.
}
\item{method}{
Monte-Carlo approximation (mc, default) or numerical approximation (uniroot).
}
}
\details{
The considered distribution has c.s.f. (1+ gamma^(-rho/gamma))^(1/rho).%%  ~~ If necessary, more details than the description above ~~
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
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~Burr }% __ONLY ONE__ keyword per line
