plotMortality <- function(param,S){
  
  require(ggplot2)
  
  
  M2 <- S$M2[param$tEnd/param$dt,,]
  
  p1 <- ggplot(data = data.frame(x = S$w, y = M2[1,]), aes(x=x,y=y)) + geom_line(size = 1)+
    scale_y_continuous(name = 'Predation mortality')+
    scale_x_log10("weight (g)", limits= c(1,max(param$wInf)))+theme_bw()
  
  
  for (i in 2:param$nSpecies){
    p1 <- p1 + geom_line(data = data.frame(x = S$w, y = M2[i,], size = 0.5, alpha = 1/2))
  }
  p1
  
  return(p1)
  
}