gradient = function(x){

n <- length(x)  
h <- 1:length(x) # Indexvalue 
g <- matrix(0,n)  # Gradient calculation 


# Take forward differences on left and right edges

g[1] <- (x[2] - x[1])/(h[2]-h[1])
g[n] <- (x[n] - x[n-1])/(h[n]-h[n-1])

if (n > 2){
hn <- h[3:n] - h[1:(n-2)];
g[2:(n-1)] <- (x[3:n]-x[1:(n-2)])/hn;
}

return(g)

}