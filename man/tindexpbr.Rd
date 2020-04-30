\name{tindexpbr}
\alias{tindexpbr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bias-reduced version of the proportionality tail index estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function computes a bias-reduced version of the proportionality estimator of the extreme value index.
}
\usage{
tindexpbr(X, k)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This estimator is defined as 1/(1+Fbar(X,e(X,1-k/n))/(k/n)*(1/(1+r))), where Fbar is the empirical c.s.f., e(X,1-k/n) the empirical expectile of level 1-k/n and 1+r=(1-mean(X)/e(X,1-k/n))*n/(n-2*k)/(1+b Fbar(X,e(X,1-k/n))^{(-rho)}/(1-rho-gamma)). Note that in this case, the assumption A(t)=b gamma t^{(-rho)}, and the associated estimators of b and rho may be found in Caeiro et al. (2005).
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
Caeiro, F., Gomes, M. I., and Pestana, D. (2005). Direct reduction of bias of the classical Hill estimator. Revstat, 3(2):113-136.
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
\code{\link{tindexp}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~hill }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~expectile }% __ONLY ONE__ keyword per line
\keyword{ ~reduced-bias }% __ONLY ONE__ keyword per line
