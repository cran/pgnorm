\name{rpgnorm}
\alias{rpgnorm}

\title{A random number generator for the \eqn{p}-generalized normal distribution}
\description{The function simulates the univariate \eqn{p}-generalized normal distribution by using one of the following methods: the \eqn{p}-generalized polar method (pgenpolar), the \eqn{p}-generalized rejecting polar method (pgenpolarrej), the Monty Python method (montypython), the Ziggurat method (ziggurat) and the method of Nardon and Pianca (nardonpianca).}
\usage{rpgnorm(n, p, mean, sigma, method)}

\arguments{
\item{n}{The natural number of random variables to be simulated.}
  \item{p}{A positive number expressing the form parameter of the distribution. The default is 2. In case of the Monty Python method and the Ziggurat method, p can be chosen from \eqn{(1, \infty) \cup \{ 0.25, 0.45, 0.5, 0.6, 0.75 \} }. }
  \item{mean}{A real number expressing the expectation of the distribution. The default is 0.}
  \item{sigma}{A positive number expressing the standard deviation of the distribution. The default is \eqn{\sigma_p=p^{1/p} \, \sqrt{\Gamma(3/p)/ \Gamma(1/p)}}, the natural standard deviation of the \eqn{p}-generalized normal distribution.}
  \item{method}{A string expressing the method to be used for the simulation ("pgenpolar", "pgenpolarrej", "montypython", "ziggurat" or "nardonpianca"). The default is "nardonpianca".}
}
\value{An \eqn{n}-dimensional, real vector.}

\references{S. Kalke and W.-D. Richter (2013)."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Volume 83. Issue 4.}
\author{Steve Kalke}


\examples{y<-rpgnorm(10000,3,method="pgenpolar")}

\keyword{distribution}

