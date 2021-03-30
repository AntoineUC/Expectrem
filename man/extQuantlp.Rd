\name{extQuantlp}
\alias{extQuantlp}

\title{
Composite extreme quantile estimation
}
\description{
This function estimates extreme quantiles by using intermediate Lp-quantiles.
}
\usage{
extQuantlp(X,tau,k,p,estim="lpindex",br=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data vector.
}
  \item{tau}{
Numeric (extreme) probability with value in [0,1].
}
  \item{k}{
A vector of number of upper order statistics.
}
  \item{p}{
The chosen Lp-quantile (must be greater than 1).  
}
  \item{estim}{
Tail index estimator used : Hill estimator ("Hill") or Lp-quantile estimator ("lpindex", default).
}
\item{br}{
Bias reduced version if TRUE. Default: FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of k estimates of the quantile.
}
\references{
Girard, S., Stupfler, G., and Usseglio-Carleve, A. (2020). Composite bias-reduced Lp-quantile estimators of extreme quantiles and expectiles, preprint.
}
\author{
Antoine Usseglio-Carleve
}
\note{

}


\seealso{

}
\examples{

}

\keyword{ ~extreme }
\keyword{ ~quantile }
