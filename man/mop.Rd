\encoding{latin1}
\name{mop}
\alias{mop}
\alias{mop.MOP}
\alias{mop.RBMOP}
\alias{mop.beta}
\alias{mop.rho}
\title{
Mean of order p statistic for the extreme value index
}
\description{
This function compute mean of order p (MOP) basic statistic for the extreme value index (EVI), which is indeed a simple generalisation of the Hill estimator.
}
\usage{
mop(x, k, p, method = c("MOP", "RBMOP"))
}

\arguments{
  \item{x}{
Data vector.
}
  \item{k}{
a vector of number of upper order statistics.
}
  \item{p}{
a vector of mean order.
}
  \item{method}{
Method used, ("MOP", default)
                  and reduced-bias MOP ("RBMOP").
}
}
\details{
Basic statistics for the EVI estimation, the MOP of \eqn{U_{ik}}, where
\eqn{U_{ik}= \frac{X_{n-i+1:n}}{X_{n-k:n}} } and \eqn{X_{i:n}} are order statistics, is
\deqn{A(k)= ( \frac{1}{k} \sum^k_{i=1} U^p_{ik} )^{1/p},}  for \eqn{p \neq 0.}

The new class of MOP EVI- estimators is
\deqn{H_p(k)= (1 - A^{-p}(k))/p,}   for \eqn{p \neq 0.}
At \code{p=0} the above  MOP estimator is equal to classical Hill estimator.

Reduced bias MOP EVI-estimators is\cr
\deqn{RBA(k)=H_p(k) (1- \frac{\beta (1-p H_p(k) )}{1-\rho-p H_p(k)} (\frac{n}{k})^\rho ).}

}
\value{
a matrix of EVI estimates, corresponds to \code{k} row and \code{p} columns. When \code{Method = "RBMOP"} shape and scale second order parameters estimates are also returned.
}

\references{
Brilhante, M.F., Gomes, M.I. and Pestana, D. (2013). A simple generalisation of the Hill estimator.
\emph{Computational Statistics and Data Analysis}, \bold{57}, 518-- 535.

Beran, J., Schell, D. and Stehlik, M. (2013). The harmonic moment tail index estimator: asymptotic distribution and robustness.  \emph{Ann Inst Stat Math},  Published Online.

Gomes, M.I., Brilhante, M.F. and Pestana, D. (2013). New reduced-bias estimators of a positive extreme value index. \emph{Submitted article}.
}

\author{
B G Manjunath \email{bgmanjunath@gmail.com}, Frederico Caeiro \email{fac@fct.unl.pt}
}

\examples{

}

\keyword{MOP}
\keyword{EVI}
\keyword{reduced-bias}
\keyword{Hill}
