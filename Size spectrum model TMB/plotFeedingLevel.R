plotFeedingLevel <- function(param,S){
  
  require(ggplot2)
  
  
  f <- S$f[param$tEnd/param$dt,,]
  fc <- param$ks/(param$alpha*param$h)
  
  p1 <- ggplot(data = data.frame(x = S$w, y = f[1,]), aes(x=x,y=y)) + geom_line(size = 1)+
    scale_y_continuous(name = 'Feeding level',limits = c(0,1), breaks = c(0.2,0.4,0.6,0.8,1))+
    scale_x_log10('weight (w)')+ geom_hline(yintercept = param$f0, linetype =2)+geom_hline(yintercept = fc, linetype = 3)+
    theme_classic()
  
  
  for (i in 2:param$nSpecies){
    p1 <- p1 + geom_line(data = data.frame(x = S$w, y = f[i,], size = 0.5, alpha = 1/2))
  }
  p1
  
  return(p1)
  
}