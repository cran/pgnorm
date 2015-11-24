\name{rpgangular}
\alias{rpgangular}

\title{A random number generator for the angular distribution}
\description{The function simulates the univariate angular distribution corresponding to the \eqn{p}-generalized uniform distribution on the \eqn{p}-generalized unit circle.}

\usage{rpgangular(n,p)}

\arguments{
\item{n}{The natural number of random variables to be simulated.}
\item{p}{A positive number expressing the form parameter of the distribution. The default is 2.}
}
\value{An \eqn{n}-dimensional, real vector.}

\references{S. Kalke and W.-D. Richter (2013)."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Volume 83. Issue 4.}
\author{Steve Kalke}


\examples{y<-rpgangular(10000,3)}

\keyword{distribution}

