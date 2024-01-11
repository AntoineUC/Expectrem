\name{ef}
\alias{ef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectiles of F distribution
%%  ~~function to do ... ~~
}
\description{
This function computes the expectiles of the F distribution.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ef(probs, df1, df2, tol=1e-08, maxiter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{
Numeric vector of probabilities with values in (0,1).%%     ~~Describe \code{x} here~~
}
\item{df1}{
Degree of freedom 1. Must be strictly positive.
}
\item{df2}{
Degree of freedom 2. Must be strictly greater than 2.
}
 \item{tol}{
Tolerance in the stopping criterion. Default: 1e-08.
}
\item{maxiter}{
Maximum number of iterations. Default: 100.
}
}
\details{
The considered distribution has p.d.f. (df1/df2)^(df1/2)*y^(df1/2-1)*(1+(df1/df2)*y)^(-(df1+df2)/2)/B(df1/2,df2/2) for y>0.%%  ~~ If necessary, more details than the description above ~~
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

\keyword{ ~expectiles }
