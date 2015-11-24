\name{rpgnorm_montypython}
\alias{rpgnorm_montypython}

\title{A random number generator for the \eqn{p}-generalized normal distribution}
\description{The function simulates the univariate, central, \eqn{p}-generalized normal distribution by using the Monty Python method.}
\usage{rpgnorm_montypython(n,p)}

\arguments{
\item{n}{The natural number of random variables to be simulated.}
  \item{p}{A positive number expressing the form parameter of the distribution. The default is 2. In case of the Monty Python method, p can be chosen from \eqn{(1, \infty) \cup \{ 0.25, 0.45, 0.5, 0.6, 0.75 \} }.}
}

\value{An \eqn{n}-dimensional, real vector.}

\references{S. Kalke and W.-D. Richter (2013)."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Volume 83. Issue 4.}
\author{Steve Kalke}


\examples{y<-rpgnorm_montypython(10000,3)}

\keyword{distribution}




