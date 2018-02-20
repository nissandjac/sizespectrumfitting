
plotGrowth  <- function (param,S){
require(ggplot2)
require(deSolve)
  
t <- param$tEnd
w <- S$w

iPlot <- round(S$nSave*t/param$tEnd);
idxRange <- floor(iPlot/2):iPlot;

tRange <- seq(0,50,by= 0.02)
hbar <- param$alpha*param$f0*param$h-param$ks


if (length(hbar) == 1){
hbar = hbar * matrix(1,param$nSpecies);
}
state <- param$w0
# Ode solver 
deriv <-function(t, state, parameters) {
   with(as.list(c(state, parameters)),{
     # rate of change
      dX <- interp1(w,g,state)
       # return the rate of change
       list(dX)
     }) # end with(as.list ...
   }

refDeriv <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    dX <- hbar*state^n*(1-(state/wInf)^(1-n)) # Hbar[i] and wInf [i] is sent in
    list(dX)
  })
}
  
  
wCalc <- matrix(NA,param$nSpecies,length(tRange))
wRef <- matrix(NA,param$nSpecies,length(tRange))
for (i in 1:param$nSpecies){
g = colMeans(S$g[idxRange,i,])

SolveG = ode(y = param$w0,times = tRange, func = deriv, parms = c(g,w))
wCalc[i,] <- SolveG[,2]

SolveGref = ode(y = param$w0, times = tRange, func = refDeriv,parms = list(w = w,wInf = param$wInf[i],n = param$n,hbar = hbar[i]))

wRef[i,] = SolveGref[,2]
}

# Adjust time for birth at w = 0:


t = tRange - 1/(hbar[1]*(param$n-1)) * param$w0^(1-param$n);
#width = seq(0.5,2,length.out = param$nSpecies)

p <- ggplot(data = data.frame(x = t, y = wCalc[1,]), aes (x = x,y=y))+geom_line()+theme_bw()+
  scale_y_log10("weight (g)", limits = c(0.01*min(param$wInf),2*max(param$wInf)))+
  scale_x_continuous("time (years)", limits = c(0.1,50))
p <- p + geom_line(data = data.frame(x = t, y = wRef[1,]), aes (x = x,y=y), linetype = 2)
for (i in 2:param$nSpecies){
p <- p + geom_line(data = data.frame(x = t, y = wCalc[i,]), aes (x = x,y=y))
p <- p + geom_line(data = data.frame(x = t, y = wRef[i,]), aes (x = x,y=y), linetype = 2)
}
p  
  
  }