\name{lpindex}
\alias{lpindex}
\title{Lp quantile based tail index estimator
}
\description{
This function computes the Lp quantile based tail index estimator.
}
\usage{
lpindex(X,k,p,br=FALSE)
}
\arguments{
  \item{X}{
Data vector.
}
  \item{k}{
A vector of number of upper order statistics.
  }
  \item{p}{
The chosen Lp-quantile (must be greater than 1).  
  }
  \item{br}{
Bias-reduced version if TRUE. Default: FALSE.
  }
}
\details{

}
\value{
A vector of k estimates of the tail index.
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
\keyword{ ~hill }
\keyword{ ~quantile }
