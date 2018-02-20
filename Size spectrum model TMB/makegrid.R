makegrid = function(param){

w <- matrix()
i <- 1;
w[1] <- param$w0[[1]];

while (w[i] < param$wMax){
w[i+1] <- w[i] + param$fGridexp*w[i]
i <- i + 1
}

if (length(w) < param$nSpecies){
param$fGridexp <- param$fGridexp/5
i <- 1
rm(w)
w <- list()
w[1] <- param$w0
while (w[i] < param$wMax){
w[i+1] <- w[i] + param$fGridexp*w[i]
i <- i + 1;
}}


#dw <- gradient(w);
k <- log(1+param$fGridexp);
dw <- k*param$w0*exp((0:(length(w)-1))*k)

return(list(w = w, dw = dw))
}
