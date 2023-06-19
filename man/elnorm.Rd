\name{elnorm}
\alias{elnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of the Lognormal distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the Lognormal distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
elnorm(probs, mu=0,sigma=1,niter=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{mu}{
Meanlog parameter. 
}
\item{sigma}{
sdlog parameter. Must be strictly positive.
}
\item{niter}{
Number of Newton-Raphson iterations. Default: 50.
}
}
\details{
The considered distribution has the stochastic representation X=exp(Y), where Y is a normal distribution with mean mu and standard deviation sigma. %%  ~~ If necessary, more details than the description above ~~
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
\keyword{ ~expectile }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~Lognormal }% __ONLY ONE__ keyword per line