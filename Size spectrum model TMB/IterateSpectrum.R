#
# Iterate the spectrum$ 
# Input: param: parameters, see baseparameters$m
#        S: optional parameter with the results of a previous run which is
#           used for initial conditions$
# Output: S: results, see the end of this script$
#
# -------------------------------------------------------------------------
  #
#  4/3-09:  Introduced separate grid for PP
#  25/1-08: Added allometrically expanding grid$
#  29/1-08: Dropped Z0 as a parameter and introduced Z0prefactor and
#           Z0exponent instead
#  23/1-08: Re-enabled starvation mortality
#  18/1-08: Added parameter wMax for the maximum grid size
#  18/1-08: Added parameter Ninit for intial spectrum
#  2/1-08: Updated the psi-rule to conform with von Bertalanffy growth
#  21/11-07: Fixed psi and switched to psi-rule reproduction
#  20/11-07: Updated bioenergetic model to have 1) std$metab; 2) 50-50
#            funneling to spawning and growth$
#  26/10-07: Added vulnerabilities for each species
#  6/9-07: Introduced parameter PPmin
#  5/9-07: Added param$bFishingEcosystem
#  15/8-07: Added possibility of random recruitment, quick and dirty
#  15/8-07: Made the background spectrum semi-chemostat
#  23/5-07: Updated fysiological model to prioritize gonads
#  9/5-07: Refixed psi to a step function
#  25/4-07: Introduced option to fix N(w0), bFixedN0
#  25/3-07: Rewrote to accept param struct
#  24/3-07: Added wPPcut option where the background spectrum is cut
#  21/3-07: Made a speedup on the calculation of phi and M2
#  5/1-06: Fixed calculation of psi
#  3/1-06: Added starvation mortality, Ms0
#  3/1-06: Added option for extreme switching
#  2/1-06: Added analytical solution for logistic background
#  21/10-05: Fixed solver to be a for-next loop instead of a
#  vector operation$
#  17/10-05: Changed reproduction to gonado-somatic index based
#  (pGonado)
#  17/10-05: If rPP<-0, then the background is fixed
#  16/8-05:  Added live background spectrum
#   3/5-05:  Added check on predation calculation (if kappaPP<-0)$
#  14/2-05:  Fixed missing kappaPP on predkernelPP
#  Added the specification of an initial spectrum: Ninitial$
#
IterateSpectrum <- function(param, S){

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
for (iTime in 1:iTimeMax){
# if param$bVerbose
# percent <- iTime/iTimeMax*100;
# if (floor(percent) > floor((iTime-1)/iTimeMax*100))
# fprintf('$',iTime/iTimeMax); ##ok<CTPCT>
# if (floor(percent/10) > floor((iTime-1)/iTimeMax*10))
# fprintf('#d\n',floor(percent));
# end
# 
# end
# end
# drawnow

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
  #  M2 <- matrix(1,nSpecies,1)%*%t(M2)}}}
#   else { # with species interaction matrix
#                 # Available food:
#                 for (j in 1:nSpecies){
#                 # Available food:
#                 #         phiprey <- $$$
#                 #           (dwPP$*wPP$* sum((thetaPP(j,:)' * ones(1,nGridPP)) $* nPP,1)) $$$
# #           * predkernelPP(1:nGrid,:)' $$$
# #           + (dw$*w$* sum((theta(j,:)'*ones(1,nGrid))$*N)) *predkernel';
# 
# phiprey <- (dwPP*wPP* thetaPP[j]*nPP) %*% predkernelPP[1:nGrid,]+
#                 (dw*w* colSums((theta[,j]*matrix(1,nGrid))*N)) %*%predkernel
# # Everything in smaller life stages is free game 
# # ix <- (w < 25);
# # phiprey[ix] <- (dwPP*wPP*nPP) * predkernelPP(ix,:)' + $$$
# # (dw(ix)$*w(ix)$* sum(N(:,ix))) *predkernel(ix,ix)';
# 
# 
# preySave[j,] <- phiprey
# # Then go for the feeding level:
# f[j,] <- 1-((IntakeMax[j,]/SearchVol[j,]/(phiprey+(IntakeMax[j,]/SearchVol[j,]))));
# }
# #f(:,:) <- 0$6;
# M2 <- matrix(0,nSpecies, nGrid)
# M2PP <- matrix(0,param$nxPP, nGridPP);
# for (i in 1:nSpecies){
# M2[i,] <- colSums(((theta[i,]*(dw*SearchVol[i,]))*N*(1-f)) %*% predkernel);
#                 #for j <- 1:nSpecies
#                 #  tmp <- theta(i,j)*dw $* (1-f(j,:)) $* SearchVol $* N(j,:);
#                 #  M2(i,:) <- M2(i,:) + tmp * predkernel;
#                 #end
#                 M2PP <- M2PP + colSums(((theta[i,]%*%(dw*SearchVol[i,]))*N*(1-f)) %*% predkernelPP)
# }
# #       for i <- 1:param$nxPP
# #         M2PP(i,:) <- sum(((thetaPP(:,i)*(dw$*SearchVol)) $* N $* (1-f)) $$$
#                           #           * predkernelPP);
# #       end
# #end
# }
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
}
# --------------------------------------------------------------
# --------------------------------------------------------------
rm(S)
S <- list()
#
# Numerical parameters:
#
S$t <- seq(param$dt,param$tEnd, by = param$dt*param$iSave) # The time steps where the
                           # solution is saved
S$nSave <- nSave;      # Number of timesteps where results are saved
                         #
                           # Grid:
                             #
S$nSpecies <- nSpecies;# No$ of species (trait classes)$
S$w <- w;              # Individual weight
S$dw <- dw;            # Difference between weight classes
                           #
                           # Species specific rates calculated directly from param:
                             #
S$wMature <- wMature;  # Weight at maturation of each species (wInf)
S$SearchVol <- SearchVol; # Volumentric searh rate  (weight)
#S$Imax <- IntakeMax;  # Maximum consumption (weight)
S$Z0 <- Z0;            # Background mortality (wInf)
S$psi <- psi;
S$theta <- theta;      # Interaction matrix (if used)
S$thetaPP <- thetaPP;  # Interaction with resource (if used)
                           #
                           # Growth:
                             #
S$g <- gSave;          # Growth (time,wInf,weight)
S$f <- fSave;          # Feeding level (time,wInf,weight)
                           #S$Eaten <- Eaten;     # Amount of food eaten (time,wInf)
S$phiprey <- phiprey; # Available food (weight)
S$prey <- preySave;
S$M2PP <- M2PPSave                           #S$e <- eSave;         # Available energy (Ee) (time,wInf,weight)
                           #
                           # Mortality:
                             #
S$M2 <- M2save;        # Predation mortality (time,wInf,weight)
S$Ms <- MsSave;        # Starvation mortality (time,wInf,weight)
S$Fin <- Fin              # Fishing mortality (wInf,weight)$
                           #
                           # Reproduction & recruitment:
                             #
S$Rp <- Rpsave;        # Egg production (time,wInf) measured in mass/time$
S$R <- Rsave;          # Recruitment flux (time,wInf) measured in numbers/time
                           #
                           # Spectra:
                             #
S$Biomass <- Biomass;  # Total biomass (time,wInf)
S$SSB <- SSBmsave      # Spawning stock biomass
S$Ntot <- NtotSave;    # Community spectrum exclusing the resource spectrum (time,weight)
S$N <- Nsave;          # Species spectra (time,wInf,weight)
                           #
                           # Resource:
                             #
S$wPP <- wPP;          # Weight in the resource spectrum
S$nPP <- nPPsave;      # Resource spectrum (time,1,weight)
S$dwPP <- dwPP;
S$xPP <- xPP; 
                         
                           
                           
                           
#                            if param$bVerbose
#                            tRun <- cputime-tStart;
#                            fprintf('\nExecution time: #2d:#05$2f$\n',$$$
#                                    [floor(tRun/60) mod(tRun,60)]);
#                            end
#     
return(S)
}
