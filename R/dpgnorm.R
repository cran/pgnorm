dpgnorm <-
function(y,p,mean,sigma)
{
# A function implemented by Steve Kalke

# Description: 
# Computes the density function of the p-generalized normal distribution 
# for the real Argument "y"

# Arguments: 
# p- a positiv constant (default: p=2)
# mean- a real constant, expressing the expectation (default: mean=0)
# sigma- a positiv constant, expressing the standard deviation (default: sigma=p^(1/p)*sqrt(gamma(3/p)/gamma(1/p)))

if(missing(p)){p<-2}

if(p<=0){stop("p has to be positive")}

sigma_p<- p^(1/p)*sqrt( gamma(3/p)/gamma(1/p) )

if(missing(mean)){mean<-0}

if(missing(sigma)){sigma<-sigma_p}

if(sigma<=0){stop("sigma has to be positive")}

#scaling and shifting the argument y
x<-sigma_p/sigma*(y-mean)

#calculating the density
density<-sigma_p/sigma*( p^(1-1/p)/2/gamma(1/p)*exp(-(abs(x)^p)/p) )

return(density)

}
