\name{rpgnorm_pgenpolarrej}
\alias{rpgnorm_pgenpolarrej}


\title{A random number generator for the \eqn{p}-generalized normal distribution}
\description{The function simulates the univariate, central, \eqn{p}-generalized normal distribution by using the \eqn{p}-generalized rejecting polar method.}

\usage{rpgnorm_pgenpolarrej(n,p)}

\arguments{
\item{n}{The natural number of random variables to be simulated.}
  \item{p}{A positive number expressing the form parameter of the distribution. The default is 2.}
}

\value{An \eqn{n}-dimensional, real vector.}

\references{S. Kalke and W.-D. Richter (2013)."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Volume 83. Issue 4.}
\author{Steve Kalke}


\examples{y<-rpgnorm_pgenpolarrej(10000,3)}

\keyword{distribution}




