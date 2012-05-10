zigsetup <-
function(p,n,tol){

# A function implemented by Steve Kalke

# Description: 
# Set's up the ziggurat for the one-sided, p-generalized normal distribution
# by returning a vector with the rightmost endpoints of the first n-1 rectangles


# Arguments: 
# p-   a positiv constant (default: p=2)
# n-   a positive integer, expressing the number of rectangles(default: n=2^8)
# tol- a positive constant, defining the accuracy of evaluation (default: tol=10^(-9))

if(missing(n)){n<-2^8}

if(round(n)!=n |n<=0){stop("invalid argument for positive integer n")}

if(missing(tol)){tol<-10^(-9)}

#defining the return vector 
x<-rep(1,n)

#finding an upper bound for the n-1 return values            
x1<-10  
while( x1*2*dpgnorm(x1,p)+1-igamma(((x1)^p)/p,1/p)>= 1/n ){x1<-10*x1} 
#lower bound for the n-1 return values
x0<-0

#nested intervals method for setting up the ziggurat
while(abs(x[1])>tol){
#rightmost endpoint of the last rectangle
x[n]<-(x0+x1)/2  
x[1]<-(-1)
#surface area of every rectangle
v<-x[n]*2*dpgnorm(x[n],p)+1-igamma((x[n]^p)/p,1/p) 

#calculation of the other rightmost endpoints
for (i in seq(n-1,1,-1)){
if( gamma(1/p)*(p^(1/p-1))*(v/x[i+1]+2*dpgnorm(x[i+1],p)) > 1)break 

if(i==1){x[1]<-x[2]-v/2/(dpgnorm(0,p)-dpgnorm(x[2],p))}
else{x[i]<-(-p*log(gamma(1/p)*(p^(1/p-1))*(v/x[i+1]+2*dpgnorm(x[i+1],p))))^(1/p) }
}
#testing if the calculated coordinates are suitable with respect to the tolerance (tol)
if(x[1]>tol){x1<-x[n]}
if(x[1]<(-1*tol)){x0<-x[n]}

}

return(x[2:n])
}
