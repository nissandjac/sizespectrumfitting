# Plot the temporal evolution of biomasses in the system (check for species dying)

plotBiomass <- function(param,S){
  
  tEnd <- param$tEnd/param$dt
  Biomass <- S$Biomass[tEnd,]
  wInf <- param$wInf
  
  dfnew <- data.frame(wInf = t(wInf), Biomass = Biomass)
  
  p <- ggplot(data = dfnew, aes(x = wInf, y = Biomass))+geom_point()+
    theme_bw()+scale_y_log10('Biomass')+scale_x_log10('WInf (years)')#+geom_smooth(method = 'lm', se = F, color = 'black')
  
  p
  return(p)
}
