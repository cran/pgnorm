rpgnorm_pgenpolarrej <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the p-generalized rejecting polar method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated


k<-floor(2/p)
p1<-2/p-k
p2<-1-p1

#to generate n p-generalized normal distributed random numbers we need to generate m=floor((n+1)/2) pairs of p-generalized normal distributed random numbers
m<-floor((n+1)/2)

V<-matrix(rep(0,2*m),nrow=m)

#defining the vector of generalized radiuses
R<-rep(0,m)

#defining the return vector
Z<-rep(0,n)

#simulation of m generalized radiuses and m random vectors which follow the uniform distribution on the p-generalized circle
for (i in 1:m){
s<-sum(log(runif(k)))
if(k!=2/p){
U<-runif(2)
while (U[1]^(1/p1)+U[2]^(1/p2)>1){U<-runif(2)}
s<-s+log(runif(1))*(U[1]^(1/p1))/ ( U[1]^(1/p1) + U[2]^(1/p2)  )
     }
#generation of the generalized radius
R[i]<-(-p*s)^(1/p)
#rejection method for the simulation of the uniform distribution on the p-generalized circle
U<-runif(2)
while(U[1]^p+U[2]^p>1){U<-runif(2)}
V[i,]<-U
   }

#generation of m pairs of random signs
S1<-sample(c(-1,1),m,replace=TRUE)
S2<-sample(c(-1,1),m,replace=TRUE)

#transforming the sample of m uniformly distributed random vectors into a sample of m generalized uniform basis vectors
U<-cbind(S1*V[,1],S2*V[,2])/((abs(V[,1])^p+abs(V[,2])^p)^(1/p))

#multiplying the generalized radius with the generalized uniform basis vector
Y<-U*R

#returning n of the 2m generated p-generalized normal distributed random numbers
Z[1:m]=Y[,1]
if (n>m){
	Z[(m+1):n]=Y[1:(n-m),2]
	}
return(Z)

}
