# baserun to simulate a multispecies system
setwd("C:/Users/Nis/Dropbox/UW/Fluctuations paper/Multispecies model")

source('load_files.R')
source('calcSSB.R')
source('YieldCalc.R')
W <- 10^seq(log10(30),log10(3000),length.out = 15) # 10 species in logspace 
param <- baseparameters(W,kappa = 0.005,h = 10)
param$F0 <- rep(0.2,param$nSpecies)
param$fishing <- "Trawl" # See "fishing.R" for definitions
#param$R0 <- rep(1,param$nSpecies)
#param$R.sd <- c(0.0,0.0)

param$wRcut <- 1
param$tEnd <- 150

S0 <- IterateSpectrum(param,S = NA)
#plotSpectrum(param,S)
plotBiomasstime(param,S0) # Make sure the simulation has gone to equilibrium
# Add fishing 


param$fishing <- "Trawl" # See "fishing.R" for definitions
SF <- IterateSpectrum(param,S0) # Add S here to start at initial conditions from before 
plotBiomasstime(param,SF) # Make sure the simulation has gone to equilibrium

# Plot Spectrum w. and wo. fishing 
# Sheldon scale
tEnd <- param$tEnd/param$dt # Last time step in Iteration
#Ntot <- S$Ntot[tEnd,]
#NFtot <- SF$Ntot[tEnd,]
# p <- ggplot(data = data.frame(x=SF$w,y = Ntot,z=NFtot),aes(x=x,y=y*x^2))+geom_line()+
#   geom_line(aes(y=z*x^2), color = 'red')+  
#   scale_x_log10('weight')+scale_y_log10('size spectrum')+theme_classic()
# p

# Model it 
param$tEnd <- 1
years <- 1960:2015
nyear <- length(years)

Fvec <- rep(0.3,nyear)
Fvec[15:35] <- 0.7
Fvec[36:nyear] <- 0.2

Fvec <- 0.6*exp((-log(floor(nyear/2)/(1:nyear))^2/(2*0.2)))+0.2

Fforage <- rep(0.3,nyear)
#Fforage[15:35] <- 0.1
#Fforage[36:nyear] <- 0.2

Biomass <- matrix(NA,nyear,param$nSpecies)
SSB <- matrix(NA,nyear,param$nSpecies)
Catch <- matrix(NA,nyear,param$nSpecies)
N <- array(NA,dim = c(nyear,param$nSpecies,length(SF$w)))  # 3 is number of species
M <- matrix(NA,nyear,length(SF$w)) # at weight == 30
Fsave <- matrix(NA,nyear,length(SF$w))

nsurvey <- 10 # Number of surveys 
surv.sd <- 0.5 # CV of survey error
nsize <- 20# Number of survey bins 
nsurvey.bin <- 10^(seq(log10(1),log10(param$wInf[1]*2+1), length.out = nsize))
nsurvey.means <- (nsurvey.bin[2:20]+nsurvey.bin[1:19])/2
size <- SF$w# 10^(seq(log10(1),log10(6000),length.out = nsize)) # Survey weight bins 
q <- 1e-5
surv.sel <- (1+(S0$w/(0.01 * W[1]))^-2)^-1 # Survey selectivity 
w <- S0$w

N.survey <- array(NA,dim = c(nsurvey,param$nSpecies,nsize-1,nyear)) # There always a bin less


for (i in 1:nyear){
  
  # Constant parameters 
  param$F0 <- rep(Fvec[i], param$nSpecies)
  param$F0[1] <- Fforage[i]
  
  if (i == 1){
    SF = S0
  }
  
  SF <- IterateSpectrum(param,S = SF)
  
  #N.calc <- SurveyCatch(sel.survey,SF,1:(param$tEnd/param$dt)) # Average over the year.. fix later
  w.calc <- SF$w
  
  Biomass[i,] <- SF$Biomass[param$tEnd/param$dt,]
  SSB[i,] <- calcSSB(param,SF, param$tEnd/param$dt) # Middle of the year
  N[i,,] <- SF$N[param$tEnd/param$dt,,]
  
  #widx <- which.min((param$wInf[1]-w.calc)^2)
  M[i,] <- colMeans(SF$M2[,1,]) # Only the first 
  Fsave[i,] <- SF$Fin[1,]
  
  Catch[i,] <- YieldCalc(param,SF)
  
  # Survey index 
  
  cfactor <- mean(exp(rnorm(n = 1e5,mean = 0, sd = surv.sd)))
  surv.sel.app <- t(matrix(rep(surv.sel,param$nSpecies),nrow = length(SF$w)))
  #plot(rowSums(N.length)*surv.sel*q, type = 'o', ylim = c(1,120), log = 'y')
  # Set up the measuring vector 
  cuts <- cut(SF$w,nsurvey.bin, labels = FALSE)
  
  
  for (k in 1:nsurvey){
    
    err <- exp(rnorm(n = nsize-1,mean = 0, sd = surv.sd))/cfactor #  Error within that survey
    
    for (j in 1:param$nSpecies){
      x <- aggregate(SF$N[param$tEnd/param$dt,j,]*q*surv.sel, list(cuts), sum) # Sum into bins 
      N.survey[k,j,,i] <-  x$x*err # Currently only saving the smallest species
    }
    
  }
  
}

w <- t(matrix(rep(SF$w,nyear),nrow = length(SF$w)))
dw <- t(matrix(rep(SF$w,nyear),nrow = length(SF$w)))

wsurv <- matrix(rep(nsurvey.means,nyear),nrow = length(nsurvey.means))

par(mfrow = c(2,3), mar = c(3,4,1,1))
plot(years,Fvec, type = 'o', ylab = 'fishing mortality', col = 'red',ylim = c(0,1))
lines(years,Fforage, type = 'o', col = 'black')


plot(years,Biomass[,1]/mean(Biomass[,1]), type = 'l', ylim = c(0.05,2.2), ylab = 'normalized biomass')
lines(years,Biomass[,2]/mean(Biomass[,2]), col = 'red')


plot(years,rowMeans(M), type = 'l', ylab = 'natural mortality') #w.calc > 1 & w.calc < W[1]

plot(years,Catch[,1]/mean(Catch[,1]), type = 'l', ylab = 'normalized catch')
lines(years,Catch[,2]/mean(Catch[,2]), col = 'red')

plot(years,rowSums(N[,1,]*w*dw)/mean(rowSums(N[,1,]*w*dw)), ylim = c(0.1,3), type = 'l',
     ylab = 'abundance-survey')
lines(years,rowSums(N[,2,]*w*dw)/mean(rowSums(N[,2,]*w*dw)), col = 'red', type = 'o')
# Add surveys
for (i in 1:nsurvey){
  lines(years,colSums(N.survey[i,1,,]*wsurv)/mean(colSums(N.survey[i,1,,]*wsurv)), col = alpha('black', 0.2))
  lines(years,colSums(N.survey[i,2,,]*wsurv)/mean(colSums(N.survey[i,2,,]*wsurv)), col = alpha('red', 0.2))
}

# plot the end point mortality 
plot(w.calc[w.calc > 1 & w.calc < W[1]],SF$M2[5,1,w.calc > 1 & w.calc < W[1]]+param$mu0prefactor*param$wInf[1]^param$mu0exponent, 
     log = 'xy', ylab = 'natural mortality', xlab = 'weight', type = 'l')

# create a list of survey + M2 + true N 
ls.forage <- list(N.survey = N.survey[,1,,], N = N[,1,], 
                  M2 = M, # M is at 32g
                  C = Catch[,1],
                  SSB = SSB[,1], 
                  w = SF$w,
                  Z = M+SF$Z0[1]+Fsave,
                  F0 = Fsave,
                  M.background = SF$Z0[1]) 
save(ls.forage, file = 
       'C:/Users/Nis/Dropbox/UW/Fluctuations paper/Multispecies model/single species model/forage_simulation.RData')


# Fishing mortality 
png(paste(getwd(),'/Figures/modeloutput.png', sep = ''),
    width = 640, height = 480)

par(mfrow = c(1,2))
plot(years,Fvec, type = 'o', ylab = 'fishing mortality', col = 'red',ylim = c(0,1))
lines(years,Fforage, type = 'o', col = 'black')

wsurv <- matrix(rep(nsurvey.means,nyear),nrow = length(nsurvey.means))

plot(years,Biomass[,1]/mean(Biomass[,1]), type = 'l', ylim = c(0.05,2.5), ylab = 'normalized biomass')
lines(years,rowMeans(Biomass[,2:param$nSpecies])/mean(Biomass[,2:param$nSpecies]), col = 'red')
legend('topleft', legend = c('forage fish', 'avg predators'), col = c('black', 'red'), lty = c(1,1))

dev.off()


png(paste(getwd(),'/Figures/observations.png', sep = ''),
    width = 640, height = 480)

par(mfrow = c(1,2))

plot(years,Catch[,1]/mean(Catch[,1]), type = 'l', ylab = 'normalized catch')
lines(years,Catch[,8]/mean(Catch[,8]), col = 'red')

plot(years,rowSums(N[,1,]*w*dw)/mean(rowSums(N[,1,]*w*dw)), ylim = c(0.1,3), type = 'l',
     ylab = 'abundance-survey')
lines(years,rowSums(N[,8,]*w*dw)/mean(rowSums(N[,8,]*w*dw)), col = 'red', type = 'o')

for (i in 1:nsurvey){
  lines(years,colSums(N.survey[i,1,,]*wsurv)/mean(colSums(N.survey[i,1,,]*wsurv)), col = alpha('black', 0.2))
  lines(years,colSums(N.survey[i,8,,]*wsurv)/mean(colSums(N.survey[i,8,,]*wsurv)), col = alpha('red', 0.2))
}


dev.off()



# Calculate Fmsy for species 1 

param$F0 <- rep(0.5,param$nSpecies)
Fvec <- seq(0,1.5, length.out = 25)
msy <- rep(NA,param$nSpecies)
param$tEnd <- 50

for (j in 1:length(Fvec)){
  param$F0[1] <- Fvec[j]
  
  SFmsy <- IterateSpectrum(param,S0)
  msy[j] <- YieldCalc(param,SFmsy)[1]
}

plot(Fvec,msy, type = 'l')
print(Fvec[which.max(msy)])











