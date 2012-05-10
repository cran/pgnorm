\name{pgnorm-package}
\alias{pgnorm-package}
\alias{pgnorm}
\docType{package}
\title{The \eqn{p}-generalized normal distribution}
\description{The pgnorm-package includes routines to evaluate (cdf,pdf) and simulate the univariate \eqn{p}-generalized normal distribution with form parameter \eqn{p}, expectation \eqn{mean} and standard deviation \eqn{\sigma}. The pdf of this distribution is given by \deqn{f(x,p,mean,\sigma)=(\sigma_p/ \sigma) \, C_p \,   \exp \left( - \left( \frac{\sigma_p}{\sigma } \right)^p \frac{\left|  x-mean \right|^p}{p} \right) ,} where \eqn{C_p=p^{1-1/p}/2/\Gamma(1/p)} and \eqn{\sigma_p^2=p^{2/p} \, \Gamma(3/p)/\Gamma(1/p)}, which becomes \deqn{f(x,p,mean,\sigma)=C_p \, \exp \left( - \frac{\left| x \right|^p}{p} \right),} if \eqn{\sigma=\sigma_p} and \eqn{mean=0}. The random number generation can be realized with one of five different simulation methods including the \eqn{p}-generalized polar method, the \eqn{p}-generalized rejecting polar method, the Monty Python method, the Ziggurat method and the method of Nardon and Pianca. Additionally to the simulation of the p-generalized normal distribution, the related \eqn{p}-generalized uniform distribution on the \eqn{p}-generalized unit circle and the corresponding angular distribution can be simulated by using the functions "rpgunif" and "rpgangular", respectively.}
\details{
\tabular{ll}{
Package: \tab pgnorm\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-05-10\cr
License: \tab GPL (version 2 or later) \cr
LazyLoad: \tab yes\cr
}
}
\author{Steve Kalke <steve.kalke@uni-rostock.de>}
\references{S. Kalke and W.-D. Richter."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Available online: 08 Nov 2011.}

\keyword{ package }
\examples{y<-rpgnorm(10,3)}
