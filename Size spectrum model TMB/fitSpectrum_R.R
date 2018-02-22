# direc<-"~/GitHub/sizespectrumfitting/Size spectrum model TMB"
# setwd(direc)
source('load_files.R')
source('calcSSB.R')
source('YieldCalc.R')
source('gradient.R')
source('getPsi.R')
library(TMB)

source('getPredkernel.R')
W <- 10^seq(log10(30),log10(30000),length.out = 15) # 10 species in logspace 
param <- baseparameters(W,kappa = 0.005,h = rep(10, length(W)))
param$F0 <- rep(0,param$nSpecies)
param$fishing <- "Trawl" # See "fishing.R" for definitions
param$tEnd <- 50


# Read the ICES data 
load('ICES_catch.Rdata') # Catch
LH <- read.csv('North-Sea-data/life_history.csv')
LH <- LH[order(LH$Winf),]

df <-list(        wlength = param$wlength,
                  nspecies = param$nSpecies,
                  wInf = param$wInf,
                  w = param$w,
                  dw = param$dw, 
                  dwPP = param$dwPP,
                  wPP = param$wPP,
                  p = param$p,
                  wPPlength = length(param$wPP),
                  kappaR  = param$kappaR,
                  kR = param$kR,
                  h = param$h,
                  ks = param$ks,
                  k = rep(param$k, param$nSpecies),
                  gamma = param$gamma,
                  Sequ = param$Ninit,
                  psi = getPsi(param),
                  predkernel = getPredkernel(param)[[1]],
                  predkernelPP = getPredkernel(param)[[2]],
                  dt = param$dt,
                  n = param$n,
                  q = param$q,
                  v = rep(1,param$nSpecies),
                  alpha = rep(param$alpha, param$nSpecies),
                  Rmax = param$Rmax,
                  nF = rep(0.25, param$nSpecies),
                  muF  = rep(10, param$nSpecies),
                  Z0 = param$mu0prefactor * param$wInf^param$mu0exponent, # Changed from Z0 to mu0,
                  eRepro = rep(param$eRepro, param$nSpecies),
                  iTimeMax = param$tEnd/param$dt,
                  rR = 4,
                  lR = 0.25,
                  Fzero = rep(0, param$nSpecies)


)


parameters = list(a = 1)#

compile("fitSizespectrum.cpp")

dyn.load(dynlib("fitSizespectrum"))
obj <-MakeADFun(df,parameters,DLL="fitSizespectrum",Type = 'Fun' ,checkParameterOrder = F)

# Plot the size spectrum in the last year
Nsave <- obj$report()$Nsave
# Total size spectrum
Ntot <- colSums(Nsave[,,df$iTimeMax])

plot(df$w,Nsave[1,,df$iTimeMax], log = 'xy', type ='l',ylim = c(1e-13,10), ylab = 'size spectrum', xlab = 'weight')
lines(df$w,Ntot, lty = 2, lwd = 2)

for (i in 2:df$nspecies){
  lines(df$w,Nsave[i,,df$iTimeMax])
  
}





