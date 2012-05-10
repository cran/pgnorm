\name{rpgunif}
\alias{rpgunif}

\title{A random number generator for the \eqn{p}-generalized uniform distribution}
\description{The function simulates the bivariate, \eqn{p}-generalized uniform distribution on the \eqn{p}-generalized unit circle.}

\usage{rpgunif(n,p)}

\arguments{
\item{n}{The natural number of random vectors to be simulated.}
  \item{p}{A positive number expressing the form parameter of the distribution. The default is 2.}
}

\value{A real \eqn{n \times 2} matrix.}

\references{S. Kalke and W.-D. Richter."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Available online: 08 Nov 2011.}
\author{Steve Kalke}


\examples{y<-rpgunif(10000,3)}

\keyword{distribution}


