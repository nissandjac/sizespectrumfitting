#
# Calculate equilibrium solutions:
  #
calcequ <- function(param,w){
#
# Basic parameters:
#
S <- list()  
param$gamma <- mean(param$gamma)
param$h <- mean(param$h)
  
S$lambda <- 2 - param$n + param$q;

alphae <- (param$gamma *sqrt(2*pi)*param$sigma*param$beta^(S$lambda-2)*exp((S$lambda-2)^2*param$sigma^2/2))
S$alphae <- alphae

S$kappa <- param$h/((1/param$f0-1)*alphae)

S$alphap <- param$f0* param$h*param$beta^(2*param$n-param$q-1) *exp((2*param$n*(param$q-1)-param$q^2+1)*param$sigma^2/2)

S$a <- S$alphap/(param$alpha*param$h*param$f0-param$ks)
#
# Realized Winf:
  #
S$wInfRealized <- param$wInf;# * param$f0^(1/(1-param$n));
#
# Prefactor for species spectra:
  #
wInf <- S$wInfRealized 

#kappa <- wInf$^(2*param$n-param$q-3+S$a) $* dwInf ;
#
# Species spectra:
  #
S$Z0 <- param$mu0prefactor * param$wInf^param$mu0exponent; # Changed from Z0prefactor to mu0prefactor, 27/2 nsj, same for exponent
if (length(param$alpha)== 1){
alph <- param$alpha * matrix(1,param$nSpecies)
}else{
  alph <- param$alpha;
}
# Update alpha to take std$ metabolism into account assuming p<-n:
hbar <- alph*param$h*param$f0-param$ks
#
# This is valid assuming logarithmically distributed Winf:
  #
if (length(S$a) == 1){
S$a <- S$a * matrix(1,param$nSpecies)
}

S$Nequ <- matrix(0,param$nSpecies,length(w))

for (iSpecies in 1:param$nSpecies){
S$Nequ[iSpecies,] <-  param$kappaPP*wInf[iSpecies]^(-S$lambda+param$n+S$a[iSpecies])*w^(-param$n - S$a[iSpecies])*
(1 - (w/wInf[iSpecies])^(1-param$n))^((S$a[iSpecies]+wInf[iSpecies]^(1-param$n)*S$Z0[iSpecies]/hbar[iSpecies])/(1-param$n));

S$Nequ[iSpecies, w >= wInf[iSpecies]] <- 0
}


# Normalization based on constant calculated in the addendum to A&B (2006):
#
const <- NA
for (iSpecies in 1:param$nSpecies){
const[iSpecies] <- S$a[iSpecies] * gamma( (param$q+4-2*param$n)/(1-param$n) ) / 
( gamma(S$a[iSpecies]/(1-param$n)+1) * gamma( (S$a[iSpecies]+2*param$n-param$q-4)/(param$n-1) ) )

S$Nequ[iSpecies,] <- S$Nequ[iSpecies,]*const[iSpecies];
}

#
# Normalization by fitting to expected community spectrum:
  #
if (param$nSpecies > 4){
const <- colSums(S$Nequ) / (param$kappaPP*w^(-2-param$q+param$n))

wtemp <- sort(wInf);
const <- mean(const[w>=wtemp[2] & w<=wtemp[length(wtemp)-1]]);
S$Nequ <- S$Nequ/const;
}

return(list = S)
}