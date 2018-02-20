plotSpectrum <- function(param,S){
  
require(ggplot2)
w <- S$w
wInf <- param$wInf

Ntot <-  S$Ntot[param$tEnd/param$dt,w < max(param$wInf)*2]
NPP <- S$nPP[param$tEnd/param$dt,,]


p1 <- ggplot(data = data.frame(x = S$w[w < max(param$wInf)*2], y = Ntot), aes(x=x,y=y)) + geom_line(size = 2)+
  scale_y_log10('Size spectrum (#/g/M^2)',expand = c(0,0), limits = c(min(Ntot[Ntot > 0]),max(NPP[S$wPP > 1e-2])))+
  scale_x_log10('weight (g)',expand = c(0,0),limits  = c(1e-2,max(param$wInf)))+theme_classic()+
  theme(text = element_text(size=20))
N <- S$N[param$tEnd/param$dt,,]

for (i in 1:param$nSpecies){
  p1 <- p1 + geom_line(data = data.frame(x = S$w[w <param$wInf[i]], 
                                         y = N[i,w<wInf[i]]), size = 0.5, alpha = 1/2)
}
p1
p1 <- p1+ geom_line(data = data.frame(x = S$wPP, y = NPP), size = 2, linetype = 2)
return(p1)

}