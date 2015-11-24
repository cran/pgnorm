\name{zigsetup}
\alias{zigsetup}

\title{A function for setting up the Ziggurat.}
\description{The function approximates the rightmost x-coordinates of the first n-1 rectangles defining the Ziggurat in case of the central, \eqn{p}-generalized normal distribution.}
\usage{zigsetup(p, n, tol)}

\arguments{
  \item{p}{A positive number expressing the form parameter of the distribution. The default is 2. In case of the Ziggurat method, p can be chosen from \eqn{(1, \infty) \cup \{ 0.25, 0.45, 0.5, 0.6, 0.75 \} }.}
  \item{n}{The number of rectangles that build up the Ziggurat. The default is \eqn{2^8}.}
  \item{tol}{A positive number expressing the approximation accuracy of the function. The default is \eqn{10^{-9}}.}
}

\value{An \eqn{(n-1)}-dimensional, real vector.}
\references{S. Kalke and W.-D. Richter (2013)."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Volume 83. Issue 4.}
\author{Steve Kalke}

\examples{ y<-zigsetup(3,20,10^(-6)) }


\keyword{distribution}

