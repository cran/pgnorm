rpgnorm_pgenpolar <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the p-generalized polar method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

k<-floor(2/p)
p1<-2/p-k
p2<-1-p1

#to generate n p-generalized normal distributed random numbers we need to generate m=floor((n+1)/2) pairs of p-generalized normal distributed random numbers
m<-floor((n+1)/2)

#defining the vector of generalized radiuses
R<-rep(0,m)

#defining the return vector
Z<-rep(0,n)


#generation of m generalized radiuses
for (i in 1:m){
s<-sum(log(runif(k)))
if(k!=2/p){
U<-runif(2)
while (U[1]^(1/p1)+U[2]^(1/p2)>1){U<-runif(2)}
s<-s+log(runif(1))*(U[1]^(1/p1))/ ( U[1]^(1/p1) + U[2]^(1/p2)  )
     }
R[i]<-(-p*s)^(1/p)
   }

#generation of m generalized uniform basis vectors
V<-rpgunif(m,p)

#multiplying the generalized radius with the generalized uniform basis vector
Y<-V*R

#returning n of the 2m generated p-generalized normal distributed random numbers
Z[1:m]=Y[,1]
if (n>m){
	Z[(m+1):n]=Y[1:(n-m),2]
	}

return(Z)
}
