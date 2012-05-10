rpgunif <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the bivariate p-generalized uniform distribution on the p-generalized circle

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random vectors to be simulated

#generation of n random angles which follow the angular distribution corresponding to the p-generalized uniform distribution

if(missing(p)){p<-2}
if(p<=0){stop("p has to be positive")}

Phi<-rpgangular(n,p)

#simulating the p-generalized uniform distribution by transforming random angles with the p-generalized trigonometric functions
U<-cbind(cos(Phi)/((abs(cos(Phi))^p+abs(sin(Phi))^p)^(1/p)),sin(Phi)/((abs(cos(Phi))^p+abs(sin(Phi))^p)^(1/p)))

return(U)
}
