getPredkernel <- function(param){
  nGrid <- length(param$w)
  nGridPP <- param$nGridPP
  
  
  predkernel <- matrix(0,nGrid, nGrid);
  predkernelPP <- matrix(0,nGrid, nGridPP);
  
  for (j in 1:nGrid){  # loop over predator lengths
    predkernel[j,] <- exp(-0.5*(log(param$w[j]/(param$beta*param$w))/param$sigma)^2);
    predkernel[j, param$w >= param$w[j]] <- 0
    
    predkernelPP[j,] <- exp(-0.5*(log(param$w[j]/(param$beta*param$wPP))/param$sigma)^2);
  }
  
return(list(predkernel = predkernel, predkernelPP = predkernelPP))  
}