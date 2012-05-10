ppgnorm <-
function(y,p,mean,sigma)
{

# A function implemented by Steve Kalke

# Description: 
# Computes the distribution function of the p-generalized normal distribution 
# for the real Argument "y" 

# Arguments: 
# p- a positiv constant (default: p=2)
# mean- a real constant, expressing the expectation (default: mean=0)
# sigma- a positiv constant, expressing the standard deviation (default: sigma= p^(1/p)*sqrt(gamma(3/p)/gamma(1/p)) )

# Note that "igamma" from Rmetrics-library is used to evaluate the incomplete 
# gamma function \Gamma(a,x)=\int\limits_x^{\infty} t^{a-1} e^{-t} dt
# igamma(x,a)= frac{1}{\Gamma(a)}  * \int_0^x e^{-t} t^{a-1} dt 


# P(X<x)=1/2 + sign(x)* \frac{1}{2}* \frac{\Gamma(1/p)-\Gamma(1/p,\vert t \vert^p /p) }{\Gamma(1/p)}

if(missing(p)){p<-2}

if(p<=0){stop("p has to be positive")}


sigma_p<- p^(1/p)*sqrt( gamma(3/p)/gamma(1/p) )


if(missing(mean)){mean<-0}

if(missing(sigma)){sigma<- sigma_p }

if(sigma<=0){stop("sigma has to be positive")}

#scaling and shifting the argument y
x<-sigma_p/sigma*(y-mean)

#calculation of the cdf
erg<-1/2+sign(x)*1/2*igamma(abs(x)^p/p,1/p)

return(erg)
}
