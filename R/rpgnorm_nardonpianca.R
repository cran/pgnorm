rpgnorm_nardonpianca <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the method of Nardon and Pianca

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

l<-floor(1/p)
p1<-1/p-l
p2<-1-p1

#defining the return vector
R<-rep(0,n)

for (i in 1:n){

#generation of n gamma distributed random numbers with the method of Jöhnk
s<-sum(log(runif(l)))
if(l!=1/p){
U<-runif(2)
while (U[1]^(1/p1)+U[2]^(1/p2)>1){U<-runif(2)}
s<-s+log(runif(1))*(U[1]^(1/p1))/ ( U[1]^(1/p1) + U[2]^(1/p2) )
     }
#transforming a gamma distributed random number into a random number following the absolute p-generalized normal distribution
R[i]<-(-p*s)^(1/p)

}

#multiplying with random signs to transform a sample of the absolute p-generalized normal distribution into a sample of the p-generalized normal distribution 
return(R*sample(c(-1,1),n,replace=TRUE))

}
