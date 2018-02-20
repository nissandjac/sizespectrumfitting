getPsi <- function(param){
  psi <- array(0, dim = c(param$nSpecies,length(param$w)))
  
  for (iSpecies in 1:param$nSpecies){
  tmp <- param$w/param$wInf[iSpecies];
  psi[iSpecies,] <- tmp^(1-param$n)*1/(1+(tmp/param$alphaMature)^(-10));     # cut off before maturation
  # param to mature cutoff
  co <- 0.1; # 0$1 originally
  psi[iSpecies, param$w<co*param$alphaMature*param$wInf[iSpecies]] <- 0
  psi[iSpecies,tmp>1] <- 1
  
  }
  
return(psi)  ##
}