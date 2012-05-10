rpgangular <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the univariate angular distribution corresponding to the (bivariate) p-generalized uniform distribution

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if(missing(p)){p<-2}
if(p<=0){stop("p has to be positive")}


V<-matrix(rep(0,2*n),nrow=n)

#rejection method for sampling from the uniform distribution of the p-generalized circle
for (i in 1:n){
U<-runif(2)
while(U[1]^p+U[2]^p>1){U<-runif(2)}
V[i,]<-U 
}

#transforming the sample of n uniformly distributed random vectors into a sample of n random angles by using n random signs and n uniformly on {0,1} distributed random numbers
S1<-sample(0:1,n,replace=TRUE)
S2<-sample(c(-1,1),n,replace=TRUE)

return(S2*atan(V[,1]/V[,2])+pi/2+pi*S1)

}
