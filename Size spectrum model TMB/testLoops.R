phiprey2 <-  t(dwPP*wPP*nPP) %*% t(predkernelPP) + 
  (dw*w*colSums(matrix(rep(param$v,nGrid),nrow = nSpecies)*N)) %*% t(predkernel)


M2[jPred,] = param$v[jPred]*colSums((((matrix(1,nSpecies)%*%dw)*SearchVol) * N 
                                     * (1-f)) %*% predkernel)


for (jPred in 1:nSpecies){
  for (i in 1:nGrid){
    M2[jPred,i] = param$v[jPred]*colSums((((matrix(1,nSpecies)%*%dw)*SearchVol) * N * (1-f)) %*% predkernel)
    
  }
}

fNtest <- array(0, dim = c(nSpecies,nGrid))

for(k in 1:nSpecies){ # Loop over species
  for(j in 1:nGrid){ #Calculate the available N)
  fNtest[k,j] = N[k,j]*(1-f[k,j])*dw[j]*SearchVol[k,j]
  }
}
