plotBiomassDens <- function(param,S){
  
  require(ggplot2)
  w <- S$w
  wInf <- param$wInf
  
  Ntot <-  S$Ntot[param$tEnd/param$dt,w < max(param$wInf)*2]*S$w[w < max(param$wInf)*2]
  
  p1 <- ggplot(data = data.frame(x = S$w[w < max(param$wInf)*2], y = Ntot), aes(x=x,y=y)) + geom_line(size = 2)+
    scale_y_log10()+scale_x_log10()+theme_bw()
  N <- S$N[param$tEnd/param$dt,,]
  
  for (i in 1:param$nSpecies){
    p1 <- p1 + geom_line(data = data.frame(x = S$w[w <param$wInf[i]], 
                                           y = N[i,w<wInf[i]]*w[w < wInf[i]]), size = 0.5, alpha = 1/2)
  }
  p1
  
  return(p1)
  
}