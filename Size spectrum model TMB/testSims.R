S <- NA
if (length(S) == 1) {
  S <- list()
}

source("makegrid.r")
source('gradient.r')
source('fishing.r')
# ---------------------------------------------------------
# Set up grid for weights (w)
# ---------------------------------------------------------
#nGrid <- param$nGrid;
#w <- logspace(log10(param$w0), log10(param$wMax), nGrid);
#dw <- gradient(w);

wl <- makegrid(param);
w <- wl[[1]]
dw <- wl[[2]]

nGrid <- length(w);

# wPP is the primary spectrum grid

nGridPP <- param$nGridPP;


#wend <- fzero(@(x) x+param$fGridexp*x^0$75-param$w0, param$w0);
if (nGridPP == 1){
  wPP <- param$wRcut;
  dwPP <- 1/wPP
}else{
  # if (length(param$wStart) == 0){
  param$wStart <- param$w0/(param$beta*4);
  wPPtemp <- exp(log(10)*seq(log10(param$wStart),log10(param$wRcut),length.out=nGridPP+1))
  dwPPtemp <- gradient(wPPtemp);
  wPP <- wPPtemp[1:nGridPP] + dwPPtemp[1:nGridPP]/2
  dwPP <- gradient(wPP)
}

# ---------------------------------------------------------
# Primary production spectrum:
# ---------------------------------------------------------
if (nGridPP == 1){
  NinfPP <- param$kappaR
}else{
  NinfPP <- param$kappaR * wPP^param$kR
}
nPP <- NinfPP; # Start with a saturated background spectrum

rrPP <- param$rR*wPP^(-param$lR) # Weight specific growth rate

# Copy the spectrum in the x-direction if needed:
xPP <- matrix(1,nGridPP)


# if (is.na(S$nPP[param$tEnd/param$dt,1,1]) == FALSE){
#  nPP <- (S$nPP[length(wPP),,]) # Fix this once the model is running
# }

#---------------------------------------------------------
# Set up the species:
# ---------------------------------------------------------
nSpecies <- param$nSpecies;
wend <- length(w)
#
# First distribute the parameters over all species and calculate
# intermediate variable:
#
onez <- matrix(1,nSpecies)

# Preallocate some stuff 

psi <- matrix(0,nSpecies,length(w))
e <- matrix(0,nSpecies,nGrid)   
gg <- matrix(0,nSpecies,nGrid)
Fin <- matrix(0,nSpecies,nGrid)
# Maximum intake:
if (length(param$h > 1)){
  IntakeMax <- matrix(0,nSpecies,wend)
  # Search volume:
  SearchVol <- matrix(0,nSpecies,wend)
}
# Bioenergetics
Activity <- param$k * w;
StdMetab <- matrix(0,nSpecies,wend)
for (i in 1:nSpecies){
  if (length(param$ks) == 1){
    StdMetab[i,] <- param$ks * w^param$p}
  else{
    StdMetab[i,] <- param$ks[i]* w^param$p}
}
# if (length(param$ks > 1)){
# 
#   StdMetab <-  t(t((matrix(1,param$nSpecies)%*%w)^param$p)%*%(diag(ks)))
#   
# }else{
#  
# }
# Fix if rho is not allocated for all species:
if (length(param$rho)<nSpecies){
  param$rho <- param$rho(1)*onez;
  param$R0 <- param$R0(1)*onez;
}

if (length(param$alpha)==1){
  param$alpha <- param$alpha*onez
}

#if (length(param$bFixedN0)<-<-1)
#  param$bFixedN0 <- param$bFixedN0*onez;
#  param$fRandomRecruitment <- param$fRandomRecruitment*onez;
#end

Winf <- param$wInf;
wMature <- Winf*param$alphaMature;
Z0 <- param$mu0prefactor * Winf^param$mu0exponent; # Changed from Z0 to mu0

if (length(param$eRepro) == 1){
  param$eRepro <- param$eRepro * matrix(1,nSpecies)
}

if (length(S) == 0){
  if (length(param$Ninit) > 0){
    N <- param$Ninit;
  }}

for (iSpecies in 1:nSpecies){
  #
  # Initial spectrum:
  #
  
  
  # 
  # if (length(S) > 0){
  #   Ntemp <- approx(S$w, S$N[param$tEnd/param$dt,iSpecies,],w) # Might be unneccesary
  #   N[iSpecies,] <- Ntemp$y
  # }#}
  
  #  Define species specific h and gamma 
  if (length(param$h) > 1){
    IntakeMax[iSpecies,] <- param$h[iSpecies]* w^param$n;
  }else{
    IntakeMax <- param$h * w^param$n;
    #IntakeMax <- t(replicate(wend,IntakeMax))
  }
  
  if (length(param$gamma) > 1) {
    SearchVol[iSpecies,] <- param$gamma[iSpecies] * w^param$q
  }else{
    SearchVol <- param$gamma * w^param$q
    #SearchVol <- t(replicate(wend,SearchVol))
  }
  
  ##
  # Allocation to reproductive mass:
  #
  
  
  tmp <- w/param$wInf[iSpecies];
  psi[iSpecies,] <- tmp^(1-param$n)*1/(1+(tmp/param$alphaMature)^(-10));     # cut off before maturation
  # param to mature cutoff
  co <- 0.1; # 0$1 originally
  psi[iSpecies, w<co*param$alphaMature*param$wInf[iSpecies]] <- 0
  psi[iSpecies,tmp>1] <- 1
  ##
  # Fishing mortality is specified either through wFstart and F0 or
  # as ecosystem fishing through F0
  #
  if (length(param$F0 == 1)){
    param$F0 <- param$F0 * matrix(1,nSpecies)
  }
  
  if (length(param$fishing) > 0){
    type = param$fishing
    Fin[iSpecies,] <- fishing(param,iSpecies,w,type)
  }
  
  
  
  if (length(param$wFstart) > 0){
    
    if (length(param$wFstart) == 1){
      param$wFstart <- param$wFstart * matrix(1,length(param$wInf));
    }
    
    Fin[iSpecies,w < param$wFstart[iSpecies]] <- 0
  }
  
  
  if (length(param$wFend) > 0){
    
    if (length(param$wFend) == 1){
      param$wFend <- param$wFend * matrix(1,length(param$wInf));
    }
    
    Fin[iSpecies,w > param$wFend[iSpecies]] <- 0
  }
}



#-----------------------------------------------------------------------------
# If there was a previous run, override some of the defaults given above 
# -----------------------------------------------------------------------------
if (length(S) > 1){ 
  nPP <- S$nPP[dim(S$nPP)[1], , ]
  N <- S$N[dim(S$N)[1], , ]
}

# Allocate the variables:
A <- matrix(0,nSpecies,nGrid);
B <- matrix(0,nSpecies,nGrid);
S <- matrix(0,nSpecies,nGrid);
R <- matrix(0,nGrid,1);
f <- matrix(0,nSpecies, nGrid);

# -----------------------------------------------------------------------------
# Set up vars for calculating predation$
# predkernel(iPredator, iPrey)
# -----------------------------------------------------------------------------
predkernel <- matrix(0,nGrid, nGrid);
predkernelPP <- matrix(0,nGrid, nGridPP);

for (j in 1:nGrid){  # loop over predator lengths
  predkernel[j,] <- exp(-0.5*(log(w[j]/(param$beta*w))/param$sigma)^2);
  predkernel[j, w >= w[j]] <- 0
  
  predkernelPP[j,] <- exp(-0.5*(log(w[j]/(param$beta*wPP))/param$sigma)^2);
} # Does not work if nGridPP == 1
#
# Set up species interaction matrix, theta:
#
if (length(param$theta) == 0){
  theta <- 1
  thetaPP <- 1
}else{
  theta <- param$theta;
  thetaPP <- param$thetaPP;
}
#
# If species have an 'x' trait use that to set up interaction matrix:
#
# if isfield(param, 'x')
# theta <- zeros(nSpecies);
# thetaPP <- zeros(nSpecies,param$nxPP);
# for j <- 1:nSpecies
# theta(j,:) <- exp( -0$5*(param$x(j)-param$x)$^2/(param$sigmax^2) );
# thetaPP(j,:) <- exp( -0$5*(param$x(j)-xPP)$^2/(param$sigmax^2) );
# end
# else
# param$nxPP <- 1;
# #thetaPP <- NaN;
# end

#---------------------------------------------------------
# Initialize:
# ---------------------------------------------------------
iTimeMax <- param$tEnd/param$dt
nSave <- floor(param$tEnd/(param$dt*param$iSave))
Nsave <- array(0,dim = c(nSave, nSpecies, nGrid))
NtotSave <- matrix(0,nSave, nGrid)
Rsave <- matrix(0,nSave,nSpecies)
M2save <- array(0,dim = c(nSave, nSpecies, nGrid))
Z <- matrix(0,nGrid,1)
Rpsave <- matrix(0,nSave,nSpecies)
SSBmsave <- array(0, dim = c(nSave, nSpecies, nGrid))
Biomass <- matrix(0,nSave,nSpecies)
Rtotsave <- matrix(0,nSave,1)
fSave <- array(0,dim =c(nSave, nSpecies, nGrid))
#eSave <- zeros(nSave, nSpecies, nGrid);
gSave <- array(0, dim = c(nSave, nSpecies, nGrid))
nPPsave <- array(0, dim=c(nSave, 1, nGridPP))

MsSave <- array(0, dim = c(nSave, nSpecies, nGrid))
dt <- param$dt

preySave <- matrix(0,param$nSpecies,nGrid)
M2PPSave <- matrix(NA,nSave,nGridPP)

f <- matrix(0,nSpecies, nGrid)
SSBm <- matrix(0,nSpecies,nGrid)
# -------------------------------------------------------------------------
# Calculate total spectrum  (just for the first iteration):
# -------------------------------------------------------------------------

Ntot <- colSums(N);
idx <- 2:nGrid
#---------------------------------------------------------
# main loop:
# ---------------------------------------------------------
iTime <- 1

  if (theta == 1){ # no species preference matrix:
    # Available food:
    phiprey <-  t(dwPP*wPP*nPP) %*% t(predkernelPP[1:nGrid,]) + (dw*w*colSums(matrix(rep(param$v,nGrid),nrow = nSpecies)*N)) %*% t(predkernel)
    
    
    
    
    # Feeding level:
    if (length(param$h) == 1){
      for (jSpecies in 1:nSpecies){
        f[jSpecies,] <- 1-
          ( t(IntakeMax/SearchVol/(phiprey+(IntakeMax/SearchVol)) ))
      }
      
      #f(:,:) <- 0$6
      
      # Predation mortality
      M2 <- matrix(0, nGrid)
      M2PP <- matrix(0, nGridPP)
      for (jPred in 1:nSpecies){
        tmp <- dw * (1-f[jPred,]) * SearchVol * N[jPred,]
        M2 <- M2 + as.numeric(matrix(tmp,1) %*% predkernel)
        #M2PP <- M2PP + as.numeric(matrix(tmp,1) %*% predkernelPP)
      }
      M2 <- matrix(1,nSpecies,1)%*%t(M2)
      M2PP <- colSums((matrix(1,nSpecies)%*%(dw*SearchVol) * N * (1-f)) %*% predkernelPP)
      
    }else{
      for (jSpecies in 1:nSpecies){
        f[jSpecies,] <- 1-
          ( t(IntakeMax[jSpecies,]/SearchVol[jSpecies,]/(phiprey+(IntakeMax[jSpecies,]/SearchVol[jSpecies,])) ))
        
        
        #f(:,:) <- 0$6
      }
      # Predation mortality
      M2 <- matrix(0,nSpecies ,nGrid)
      M2PP <- matrix(0, nGridPP)
      for (jPred in 1:nSpecies){
        M2[jPred,] = param$v[jPred]*colSums((((matrix(1,nSpecies)%*%dw)*SearchVol) * N * (1-f)) %*% predkernel)
        
        # tmp <- dw * (1-f[jPred,]) * SearchVol[jPred,] * N[jPred,]
        
      }
      M2PP <- colSums(((matrix(1,nSpecies)%*%dw)*SearchVol * N * (1-f)) %*% predkernelPP)
    }
  }
 
  # -----------------------------------------------------------------
  # Iterate each species:
  # -----------------------------------------------------------------
  for (i in 1:nSpecies){
    #
    
    # Calc$ assimilated intake:
    #
    if (length(param$h) ==1){
      e[i,] <- param$alpha[i]*f[i,]*IntakeMax
    }else{
      e[i,] <- param$alpha[i]*f[i,]*IntakeMax[i,]
    }
    #
    # Subtract standard metabolism and activity:
    
    
    e[i,] <- e[i,] - StdMetab[i,] - Activity
    
    #
    # Starvation mortality is imposed if the available energy e<0:
    #
    Ms <- matrix(0,nGrid,1) #Add if adding starvation mortality 
    ix <- e[i,]<0;
    Ms[ix] <- -param$muS0*e[i,ix]/w[ix]
    e[i,ix] <- 0
    #
    # Calculate the energy used for spawning:
    #
    SSBm[i,] <- psi[i,]*e[i,] # Psi-rule for spawning:
    #
    # $$$and use the rest for somatic growth:
    #
    gg[i,] <- e[i,] - SSBm[i,]
    #
    # Total mortality:
    # 
    Z <- Z0[i] + M2[i,] + Ms + Fin[i,]
    #
    # Set up matrix for derivative:
    #
    A[i,idx] <- -gg[i,(idx-1)]*dt/dw[idx]
    B[i,idx] <- 1 + gg[i,idx]*dt/dw[idx] + Z[idx]*dt
    S[i,idx] <- N[i,idx]
    #
    # BC at upstream end (recruitment)
    #
    Egg <- sum(param$eRepro[i]*SSBm[i,]*N[i,]*dw)   # Egg production (mass/time)
    Rp <- 0.5*param$rho[i]*Egg/w[1]# * param$R0[i]*exp(rnorm(n = 1, mean = 0, sd = param$R.sd[i])) # Egg flux (numbers/time)
    
    # switch param$nRecruitmentType Maybe fix if needed
    # case 0,     # Physiological
    # R <- Rp;    # recruitment <- egg production
    # B(i,1) <- 1 + gg(i,1)$*dt$/dw(1) + Z(1)*dt;	
    # N(i,1) <- (N(i,1)+R*dt/dw(1))/B(i,1);            
    # case 1,     # Fixed
    # R <- param$fN0(i)*(gg(i,1)*w(1));
    # N(i,1) <- param$fN0(i);
    
    # Beverton-Holt
    R <- param$Rmax[i]*Rp/(param$Rmax[i]+Rp)
    B[i,1] <- 1 + gg[i,1]*dt/dw[1] + Z[1]*dt
    N[i,1] <- (N[i,1]+R*dt/dw[1])/B[i,1]
    
    #
    # Invert matrix
    #
    for (j in 2:nGrid){
      N[i,j] <- (S[i,j]-A[i,j]*N[i,j-1])/B[i,j]
    }
    
    #
    # save results:
    #
    
    if ((iTime %% param$iSave) == 0) {
      iSave <- iTime/param$iSave;
      #eSave(iSave,i,:) <- e(i,:);
      gSave[iSave,i,] <- gg[i,]
      Nsave[iSave,i,] <- N[i,]
      Rsave[iSave,i] <- R; # recruitment
      Biomass[iSave,i] <- sum(N[i,]*w*dw);
      Rpsave[iSave,i] <- Rp; # egg production
      SSBmsave[iSave,i,] <- SSBm[i,]
      MsSave[iSave,i,] <- Ms;
      Rtotsave[iSave] <- Rtotsave[iSave] + R;
      fSave[iSave, , ] <- f
      M2PPSave[iSave,] <- M2PP
      #Eaten[iSave,i] <- sum(dw*f[iSpecies,]*IntakeMax[iSpecies,]*N[iSpecies,]);
    }
  }
  #
  # Update the background spectrum:
  #
  
  # for (i in 1:param$nxPP){
  tmp <- (rrPP*NinfPP / (rrPP + M2PP))
  nPP <- tmp - (tmp - nPP)*exp(-(rrPP+M2PP)*dt)
  # }
  nPP[nPP<1e-300] <- 1e-300 # Fix points where the background has vanished
  #
  # Calculate total spectrum:
  #
  Ntot <- colSums(N)
  #
  # save results:
  #
  if ((iTime %% param$iSave) == 0){ 
    iSave <- iTime/param$iSave;
    NtotSave[iSave,] <- Ntot
    M2save[iSave, , ] <- M2
    M2PPSave[iSave,] <- M2PP
    nPPsave[iSave, , ] <- nPP
  }