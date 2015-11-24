\name{dpgnorm}
\alias{dpgnorm}
\title{A function to evaluate the \eqn{p}-generalized normal density}
\description{The function evaluates the density \eqn{f(x,p,mean,sigma)} of the univariate \eqn{p}-generalized normal distribution according to \deqn{f(x,p,mean,\sigma)=(\sigma_p/ \sigma) \, C_p \,   \exp \left( - \left( \frac{\sigma_p}{\sigma } \right)^p \frac{\left|  x-mean \right|^p}{p} \right) ,} where \eqn{C_p=p^{1-1/p}/2/\Gamma(1/p)} and \eqn{\sigma_p^2=p^{2/p} \, \Gamma(3/p)/ \Gamma(1/p)}.}
\usage{dpgnorm(y, p, mean, sigma)}
\arguments{
\item{y}{The real argument of the function.}
\item{p}{A positive number expressing the form parameter of the distribution. The default is 2.}
\item{mean}{A real number expressing the expectation of the distribution. The default is 0.}
\item{sigma}{A positive number expressing the standard deviation of the distribution. The default is \eqn{\sigma_p}.}
}
\value{A real number.}
\references{S. Kalke and W.-D. Richter (2013)."Simulation of the p-generalized Gaussian distribution." Journal of Statistical Computation and Simulation. Volume 83. Issue 4.}
\author{Steve Kalke}

\examples{y<-dpgnorm(0,3,1,2)}

\keyword{distribution}
