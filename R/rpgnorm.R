rpgnorm <-
function(n,p,mean,sigma,method){

# A function implemented by Steve Kalke

# Description: 
# Samples from the univariate p-generalized Gaussian distribution 

# Arguments: 
# n- a natural number expressing the number of random numbers to be generated
# p- a positiv constant (default: p=2)
# mean- a real constant expressing the expectation (default: mean=0)
# sigma- a positiv constant expressing the standard deviation (default: sigma= p^(1/p)*sqrt(gamma(3/p)/gamma(1/p)) )
# method- a string expressing the method to be used for the simulation ( either "nardonpianca", "montypython", "pgenpolar", "pgenpolarrej" or "ziggurat", the default is "nardonpianca")

if(missing(p)){p<-2}

if(p<=0){stop("p has to be positive")}

sigma_p<- p^(1/p)*sqrt( gamma(3/p)/gamma(1/p) )

if(missing(mean)){mean<-0}

if(missing(sigma)){sigma<- sigma_p}

if(sigma<=0){stop("sigma has to be positive")}

if(missing(method)){method<-"nardonpianca"} 

if(method !="nardonpianca" && method != "montypython" && method != "pgenpolar" && method != "pgenpolarrej" && method != "ziggurat"){stop("improper simulation method")}



#sampling from the standardized p-generalized normal distribution with respect to the chosen method

if(method=="nardonpianca"){y<-rpgnorm_nardonpianca(n,p)}
if(method=="montypython"){y<-rpgnorm_montypython(n,p)}
if(method=="pgenpolar"){y<-rpgnorm_pgenpolar(n,p)}
if(method=="pgenpolarrej"){y<-rpgnorm_pgenpolarrej(n,p)}
if(method=="ziggurat"){y<-rpgnorm_ziggurat(n,p)}

#scaling and shifting of the standardized p-generalized normal distributed sample
return(sigma/sigma_p*y+mean)
}
