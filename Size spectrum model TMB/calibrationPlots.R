# Calibration plots - Check these for succesful calibration of an LME 

calibrationPlots<- function(state_new,param,SF,Fmsy) {# State_new data, parameters from size based model, and equilibrium simulation

require(ggplot2)

# First plot the biomasses 

tEnd <- param$tEnd/param$dt
SSB <- calcSSB(param,SF,tEnd)
Biomass <- SF$Biomass[tEnd,]

SSBio <- matrix(NA,param$nSpecies)
SSBio[state_new$totFlag == 0] <- Biomass[state_new$totFlag == 0]
SSBio[is.na(SSBio)] <- SSB[which(is.na(SSBio) == 1)]

df <- data.frame(SSBIO = SSBio, bio_obs = state_new$SSBio, wInf = state_new$wInf,residuals = state_new$SSBio-SSBio)

p1 <- ggplot(data = df, aes(x = wInf,y = bio_obs)) + geom_point(size = 4)+
  geom_point(aes(y = SSBIO), color = 'blue')+
  scale_x_log10('Asymptotic weight (g)', breaks = 10^seq(log10(100),log10(1e5),length.out = 4))+
  scale_y_log10('Biomass (MT)', breaks = 10^seq(log10(1000),log10(1e6),length.out = 4))+
  theme_classic()+theme(text = element_text(size=20))

# Observed minus expected biomasses - do they allign? 
p2 <- ggplot(data = df, aes(x = log10(SSBIO),y = log10(bio_obs))) + geom_point(size = 4)+
  theme_classic()+
  scale_y_continuous(name = 'log10(Observed Biomass)')+
  scale_x_continuous(name = 'log10(Predicted Biomass)')+
  geom_abline(intercept = 0, slope = 1)

ml1 <- lm(log10(bio_obs) ~ log10(SSBio), data = df)
sm <- summary(ml1)
R_2 <- round(sm$r.squared, digits = 2)

p2 <- p2 + annotate('text', label =  R_2, 
                    y=max(log10(df$bio_obs)),x = min(log10(df$SSBIO)),hjust = -1, parse = T)

# Plot the size spectrum 

p3<- plotSpectrum(param,SF)

# Plot growth
  
p4 <- plotAvgVonB(param,SF,state_new)

# Plot the catch 

for (i in 1:param$nSpecies){
  if(state_new$Catch[i] == 0)
    state_new$Catch[i] <- NA
}
  
  
state_new$CatchSS <- state_new$Catch
state_new$CatchFlag <- NA
state_new$CatchFlag[is.na(state_new$Catch)] <- 0
state_new$CatchFlag[is.na(state_new$Catch)== 0] <- 1
state_new$CatchSS[is.na(state_new$Catch)] <- state_new$Landings[is.na(state_new$Catch)]




df <- data.frame(wInf = t(param$wInf), YieldSS = log10(YieldCalc(param,SF)), Yield_obs = log10(state_new$CatchSS),
                 grps = state_new$CatchFlag) 
df[is.na(df[,2]),] <- NA

p5 <- ggplot(df,aes(x = YieldSS,y = Yield_obs)) +geom_point(aes(group = grps))+
  theme_classic()+
  scale_y_continuous(name = 'log10(Observed Yield)')+
  scale_x_continuous(name = 'log10(Predicted Yield)')+
  geom_abline(intercept = 0, slope = 1)
# Plot mortality 
M2 <- SF$M2[tEnd,,]
M1 <- NA 
w <- SF$w

idx <- NA
for (i in 1:param$nSpecies){
  idx[i] <- which.min((param$alphaMature*param$wInf[i] - w)^2)
  
}


for (i in 1:param$nSpecies){
  # Remove the mortalities that are unused (> wInf)
  M2[i,w >param$wInf[i]] <- NA
  M1[i] <- M2[i,idx[i]]+SF$Z0[i]
}


M2plot <- ggplot(data = data.frame(x = SF$w,y = M2[1,]),aes(x = x, y = y))+geom_line()+
  scale_y_continuous(name = 'natural mortality',limits = c(0,max(M2[,w>1])))+
  scale_x_log10(name = 'weight (g)',
                limits = c(1,max(param$wInf*2)),
                breaks = c(1,10,100,1000,1e5))+
  theme_classic()

for (i in 2:param$nSpecies){
  M2plot <- M2plot + geom_line(data = data.frame(x = SF$w,y = M2[i,]),aes(x = x, y = y))
  
}

if (length(which(is.na(state_new$M) == 0)) > 0){
M2plot <- M2plot + geom_point(data = state_new, aes(x = wInf, y = M), shape = 2)
}


M1_df <- data.frame(x = t(param$wInf),y = M1)
M2plot <- M2plot + geom_point(data = M1_df, aes(x = x, y =y), shape = 3, color = 'red')


# Plot size specific exploitation 
# Plot an exploitation index averaged over different asymptotic sizes 

state_new$CatchSS <- state_new$Catch
state_new$CatchFlag <- NA
state_new$CatchFlag[is.na(state_new$Catch)] <- 0
state_new$CatchFlag[is.na(state_new$Catch)== 0] <- 1
state_new$CatchSS[is.na(state_new$Catch)] <- state_new$Landings[is.na(state_new$Catch)]

state_new$eidx <- state_new$CatchSS/state_new$SSBio

# Bin the wInf's in log groups 
cuts <- cut(state_new$wInf, 10^seq(log10(10),log10(100000),length.out = 4), 
            labels = F)

cLabels <- c('Small','Medium','Large')

dfSmall <- state_new[cuts == 1,]
dfMedium <- state_new[cuts == 2,]
dfLarge <- state_new[cuts == 3,]

dfMean <- data.frame(eidx = c(mean(dfSmall$eidx,na.rm=T),mean(dfMedium$eidx,na.rm=T),mean(dfLarge$eidx,na.rm=T)), grps = cLabels)
dfMean$grps <- reorder(dfMean$grps, c(1,2,3))

# Make a new datarame that bins the data 
p7 <- ggplot(dfMean, aes(x=as.factor(grps),y=eidx))+geom_bar(stat='identity')+theme_classic()+
  scale_x_discrete('Size group')+scale_y_continuous('Exploitation',expand = c(0,0), limits = c(0,1))


p8 <- plotRpR(param,SF)

# Plot Rmax wInf 
param$Rmax <- as.numeric(as.matrix(param$Rmax))
df <- data.frame(wInf = t(param$wInf), Rmax = param$Rmax)
p9 <- ggplot(data=df,  aes(x=wInf,y = Rmax))+scale_y_log10('Maximum recruitment')+
  scale_x_log10('asymptotic weight (g)')+theme_classic()+geom_point()

p10 <- plotBiomasstime(param,SF)
p11 <- plotFeedingLevel(param,SF)

p12 <- ggplot(data=state,aes(x=wInf,y=Fmsy))+geom_point()+
  geom_point(data=data.frame(x=t(param$wInf),y=Fmsy),aes(x=x,y=y), color='red')+
  theme_classic()+
  scale_x_log10()+scale_y_continuous(limits = c(0,1))
p12


grid.arrange(p1,p2,p3,p4,p5,M2plot,p7,p8,p9,p10,p11,p12)
}