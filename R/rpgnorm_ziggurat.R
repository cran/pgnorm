rpgnorm_ziggurat <-
function(n,p,x){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the Ziggurat method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated
# x- a real vector containing the rightmost endpoints of the 2^8 ziggurat-rectangles

#calculation of the ziggurat coordinates
if(missing(x)){x<-zigsetup(p,2^8)}

#loading the tail algorithm constants psi and beta to define a majorant of the p-generalized normal pdf in the case 0<p<1
if(p<1){       
data(datasetpgnzig)
B<-datasetpgnzig
psi<-B[,2]%*%(B[,1]==p)
beta<-B[,3]%*%(B[,1]==p)
  }

R<-rep(0,n)

#sampling n times from the absolute value of a standardized p-generalized normal distributed random variate with the ziggurat method
for (i in 1:n){
while (R[i]==0){
#choosing a rectangle uniformly from the ziggurat 
v<-sample(1:(2^8),1)
#the chosen rectangle is the last rectangle
if (v==2^8){b<-x[v-1]+(1-igamma(x[v-1]^p/p,1/p))/2/dpgnorm(x[v-1],p) 
#generating the x-coordinate of the random point from the last rectangle
u_1<-b*runif(1)
#testing if the random rectangle point is situated between the positive x-axis and the graph of the p-generalized normal pdf
if(u_1<=x[v-1]){R[i]<-u_1}
#tail algorithm 
else{

U<-runif(2)
if(p>=1){
while(U[2]>=(U[1]^(-1))*exp(x[v-1]^p/p-1/p*(x[v-1]-log(U[1])/(x[v-1]^(p-1)))^p)){U<-runif(2)}
R[i]<-x[v-1]-log(U[1])/(x[v-1]^(p-1))
   }
else{
while( U[2]*exp(-x[v-1]^p/p)*(U[1]^(beta/(beta-1)))>exp( -1/p*( ( x[v-1]+1/psi*( U[1]^(1/(1-beta))-1 ))^p) )  ){U<-runif(2)}
R[i]<-x[v-1]+1/psi*( U[1]^(1/(1-beta))-1 )
     }

}

}
#the chosen rectangle is not the last rectangle
else{
#generating the x-coordinate of the random point from the chosen rectangle
u_1<-x[v]*runif(1)

#testing if the random rectangle point is situated between the positive x-axis and the graph of the p-generalized normal pdf
if(v>1){if (u_1<=x[v-1]){R[i]<-u_1}
else{
#generating the y-coordinate of the random point from the chosen rectangle
u_2<-2*(dpgnorm(x[v-1],p)-dpgnorm(x[v],p))*runif(1)+2*dpgnorm(x[v],p)

if (u_2<=2*dpgnorm(u_1,p)){R[i]<-u_1}

}
   }
else{ 
#generating the y-coordinate of the random point from the chosen rectangle
u_2<-2*(dpgnorm(0,p)-dpgnorm(x[1],p))*runif(1)+2*dpgnorm(x[1],p)
if(u_2<=2*dpgnorm(u_1,p)){R[i]<-u_1}
     }


}
}
}

#multiplying with random signs to transform a sample of the absolute p-generalized normal distribution into a sample of the p-generalized normal distribution
R2<-sample(c(-1,1),n,replace=TRUE)*R

return(R2)

}
