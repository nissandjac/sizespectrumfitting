// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
 // Integer data
 // DATA_INTEGER(iTimeMax);
  DATA_INTEGER(wlength); // length of w
  DATA_INTEGER(wPPlength);
 // DATA_INTEGER(nlength); // Number of weight classes in survey
  DATA_INTEGER(nspecies); // Number of species

  // Vector data
  DATA_VECTOR(wInf); // Asymptotic weights
  DATA_VECTOR(w); // weight vector
  DATA_VECTOR(dw); // dw vector
  DATA_VECTOR(wPP); // Background spectrum
  DATA_VECTOR(dwPP); // Background spectrum dw
  DATA_VECTOR(h); // Growth rates
  DATA_VECTOR(ks); // Some metabolism
  DATA_VECTOR(k); // Some other metabolism
  DATA_VECTOR(gamma); // Search vol
  DATA_VECTOR(v); // Vulnerability
  DATA_VECTOR(alpha); // Energy stuff
  DATA_VECTOR(nF); // Fishing selectivity
  DATA_VECTOR(muF); // Fishing selectivity
  DATA_VECTOR(Z0); // Background mortality

//  // Array data
  DATA_ARRAY(Sequ); // Initial size distribution
  DATA_ARRAY(psi); // Mortality ogive
  DATA_MATRIX(predkernel); // Predation kernels per species
  DATA_MATRIX(predkernelPP); // Predation kernel for PP
//
//  // Scalar data
//  // DATA_SCALAR(Zbase);
//  //DATA_SCALAR(dt);
  DATA_SCALAR(n);
  DATA_SCALAR(p);
  DATA_SCALAR(q); // Search volume exponent
 // DATA_SCALAR(Rmax); // Carrying capacity
  DATA_SCALAR(kR); // Growth rate of PP
  DATA_SCALAR(kappaR); // carrying capacity of PP
  PARAMETER(a);

//
//  // Random effects (F and M)

// Calculate the species standard metabolism
  array<Type>activity(nspecies,wlength);
  array<Type>IntakeMax(nspecies,wlength);
  array<Type>SearchVol(nspecies,wlength);

array<Type>stdmetab(nspecies,wlength);

  for(int j=0;j<nspecies;j++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Loop over sizes
            stdmetab(j,i) = ks(j)*pow(w(i),p);
            activity(j,i) = k(j)*w(i);
            IntakeMax(j,i) = h(j)* pow(w(i),n);
            SearchVol(j,i) = gamma(j) * pow(w(i),q);


        }
    }

REPORT(stdmetab)
REPORT(wlength)
REPORT(activity)
REPORT(IntakeMax)
REPORT(SearchVol)


//for(int i=0;i<nspecies;i++){ // Loop over species
//    for(int j=0;j<iTimeMax;j++){ // Loop over years
//      logF(i,j)=U(i,j);
//    }
//}


vector<Type> Fin(nspecies,wlength); // Survey selectivity
//
for(int k=0;k<nspecies;k++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Initialize vectors
            // Fishing selectivity
        Fin(k,i) =pow( (Type(1.0)+pow( (w(i)/(nF(k) * wInf(k))),-muF(k))),Type(-1));
    }
}
//
//
//vector<Type> surveysel(wlength); //  Survey selectivity
//
////Type muFs = 1;
//for(int i=0;i<wlength;i++){ // Initialize vectors
//  // survey selectivity
// surveysel(i) =pow( (Type(1.0)+pow( (w(i)/(nFs * wInf)),-muFs)),Type(-1));
//
//}


//for(int i=1;i<iTimeMax;i++){ // Calculate the total mortality
//
//    for(int j=0;j < wlength;j++){
//      Fin(j,i) = Fin(j,i-1)*exp(logF(i)); // Fishing
//    }
//}


vector<Type> nPP(wPPlength); // Background spectrum
 for(int i=0;i<wPPlength;i++){ // Loop over species

  nPP(i) = kappaR * pow(wPP(i),kR);
}

array<Type> N(nspecies,wlength); // Size spectrum
for(int k=0;k<nspecies;k++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Initialize vectors
                N(k,i) = Sequ(k,i);
                }
            }

//for(int i=0;i<iTimeMax;i++){ // start time loop

    // Calculate the available food
    vector<Type>phiPP(wPPlength);
    for(int j=0;j<wPPlength;j++){ // Calculate the available PP
    phiPP(j) = dwPP(j)*wPP(j)*nPP(j);
    }
    vector<Type>Ntot(wlength);
    for(int j=0;j<wlength;j++){ // Calculate the available N
              vector<Type> m1_col1 = N.col(j);    // 1st row of m1
              Ntot(j) = m1_col1.sum();
//
    }

    vector<Type>phiN = Ntot*dw*w;

//
    vector<Type>phiprey(wlength); // Available prey
    phiprey = predkernelPP*phiPP+predkernel*phiN;

    array<Type>f(nspecies,wlength); // Feeding level

    // Calculate the feeding level

    for(int k=0;k<nspecies;k++){ // Loop over species
            for(int j=0;j<wlength;j++){ // Calculate the available N
            f(k,j) = Type(1.0)-(IntakeMax(k,j)/SearchVol(k,j)/(phiprey(j)+(IntakeMax(k,j)/SearchVol(k,j))));
            }
    }

    // Calculate mortality
    array<Type> Mtwo(nspecies ,wlength);
    vector<Type>MtwoPP(wPPlength);
    matrix<Type>fN(nspecies,wlength);

    for(int k=0;k<nspecies;k++){ // Loop over species
      for(int j=0;j<wlength;j++){ // Calculate the available N)
        fN(k,j) = N(k,j)*(1.0-f(k,j))*dw(j)*SearchVol(k,j);
      }
    }

    matrix<Type> M2tmp = fN*predkernel;
    matrix<Type> MPP = fN*predkernelPP;
    array<Type> M2(nspecies,wlength);
    array<Type>energy(nspecies,wlength);
    array<Type>SSBm(nspecies,wlength);
    array<Type>gg(nspecies,wlength);
    array<Type>Z(nspecies,wlength);

    // Population matrices for Von Voerster
    matrix<Type>A(nspecies,wlength);
    matrix<Type>B(nspecies,wlength);
    matrix<Type>S(nspecies,wlength);


    for(int k=0;k<nspecies;k++){ // Loop over species


        for(int j=0;j<wlength;j++){ // Loop over lengths
              vector<Type> m1_row1 = M2tmp.col(j);    // 1st row of m1
              M2(k,j) = v(k)*m1_row1.sum(); // Sum all contributions to mortality
              energy(k,j) = alpha(k)*f(k,j)*IntakeMax(k,j)-stdmetab(k,j)-activity(k,j); // Available energy
              SSBm(k,j) = psi(k,j)*energy(k,j); // Spawning rule
              gg(k,j) = energy(k,j)-SSBm(k,j); // Remaining energy used for growth
              Z(k,j) = Z0(k) + M2(k,j)+ Fin(k,j);
        }
// Set up matrices
    for(int idx=1;idx<wlength;idx++){ // Calculate the available N
    A(k,idx) = -Type(1.0)*gg(k,(idx-1))*dt/dw(idx);
    B(k,idx) = 1 + gg(k,idx)*dt/dw(idx) + Z(idx)*dt;
    S(k,idx) = N(k,idx);
    }

    A(k,0) = 0;
    B(k,0) = 0;
    S(k,0) = 0;

REPORT(A)
REPORT(B)
REPORT(S)

    }
//Egg <- sum(param$eRepro[i]*SSBm*N[i,]*dw)   # Egg production (mass/time)
//Rp <- 0.5*param$rho[i]*Egg/w[1]# * param$R0[i]*exp(rnorm(n = 1, mean = 0, sd = param$R.sd[i])) # Egg flux (numbers/time)
//
//# switch param$nRecruitmentType Maybe fix if needed
//# case 0,     # Physiological
//# R <- Rp;    # recruitment <- egg production// Create a file to run the TMB version of sprat
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
 // Integer data
 // DATA_INTEGER(iTimeMax);
  DATA_INTEGER(wlength); // length of w
  DATA_INTEGER(wPPlength);
 // DATA_INTEGER(nlength); // Number of weight classes in survey
  DATA_INTEGER(nspecies); // Number of species

  // Vector data
  DATA_VECTOR(wInf); // Asymptotic weights
  DATA_VECTOR(w); // weight vector
  DATA_VECTOR(dw); // dw vector
  DATA_VECTOR(wPP); // Background spectrum
  DATA_VECTOR(dwPP); // Background spectrum dw
  DATA_VECTOR(h); // Growth rates
  DATA_VECTOR(ks); // Some metabolism
  DATA_VECTOR(k); // Some other metabolism
  DATA_VECTOR(gamma); // Search vol
  DATA_VECTOR(v); // Vulnerability
  DATA_VECTOR(alpha); // Energy stuff
  DATA_VECTOR(nF); // Fishing selectivity
  DATA_VECTOR(muF); // Fishing selectivity
  DATA_VECTOR(Z0); // Background mortality

//  // Array data
  DATA_ARRAY(Sequ); // Initial size distribution
  DATA_ARRAY(psi); // Mortality ogive
  DATA_MATRIX(predkernel); // Predation kernels per species
  DATA_MATRIX(predkernelPP); // Predation kernel for PP
//
//  // Scalar data
//  // DATA_SCALAR(Zbase);
//  //DATA_SCALAR(dt);
  DATA_SCALAR(n);
  DATA_SCALAR(p);
  DATA_SCALAR(q); // Search volume exponent
 // DATA_SCALAR(Rmax); // Carrying capacity
  DATA_SCALAR(kR); // Growth rate of PP
  DATA_SCALAR(kappaR); // carrying capacity of PP
  PARAMETER(a);

//
//  // Random effects (F and M)

// Calculate the species standard metabolism
  array<Type>activity(nspecies,wlength);
  array<Type>IntakeMax(nspecies,wlength);
  array<Type>SearchVol(nspecies,wlength);

array<Type>stdmetab(nspecies,wlength);

  for(int j=0;j<nspecies;j++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Loop over sizes
            stdmetab(j,i) = ks(j)*pow(w(i),p);
            activity(j,i) = k(j)*w(i);
            IntakeMax(j,i) = h(j)* pow(w(i),n);
            SearchVol(j,i) = gamma(j) * pow(w(i),q);


        }
    }

REPORT(stdmetab)
REPORT(wlength)
REPORT(activity)
REPORT(IntakeMax)
REPORT(SearchVol)


//for(int i=0;i<nspecies;i++){ // Loop over species
//    for(int j=0;j<iTimeMax;j++){ // Loop over years
//      logF(i,j)=U(i,j);
//    }
//}


vector<Type> Fin(nspecies,wlength); // Survey selectivity
//
for(int k=0;k<nspecies;k++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Initialize vectors
            // Fishing selectivity
        Fin(k,i) =pow( (Type(1.0)+pow( (w(i)/(nF(k) * wInf(k))),-muF(k))),Type(-1));
    }
}
//
//
//vector<Type> surveysel(wlength); //  Survey selectivity
//
////Type muFs = 1;
//for(int i=0;i<wlength;i++){ // Initialize vectors
//  // survey selectivity
// surveysel(i) =pow( (Type(1.0)+pow( (w(i)/(nFs * wInf)),-muFs)),Type(-1));
//
//}


//for(int i=1;i<iTimeMax;i++){ // Calculate the total mortality
//
//    for(int j=0;j < wlength;j++){
//      Fin(j,i) = Fin(j,i-1)*exp(logF(i)); // Fishing
//    }
//}


vector<Type> nPP(wPPlength); // Background spectrum
 for(int i=0;i<wPPlength;i++){ // Loop over species

  nPP(i) = kappaR * pow(wPP(i),kR);
}

array<Type> N(nspecies,wlength); // Size spectrum
for(int k=0;k<nspecies;k++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Initialize vectors
                N(k,i) = Sequ(k,i);
                }
            }

//for(int i=0;i<iTimeMax;i++){ // start time loop

    // Calculate the available food
    vector<Type>phiPP(wPPlength);
    for(int j=0;j<wPPlength;j++){ // Calculate the available PP
    phiPP(j) = dwPP(j)*wPP(j)*nPP(j);
    }
    vector<Type>Ntot(wlength);
    for(int j=0;j<wlength;j++){ // Calculate the available N
              vector<Type> m1_col1 = N.col(j);    // 1st row of m1
              Ntot(j) = m1_col1.sum();
//
    }

    vector<Type>phiN = Ntot*dw*w;

//
    vector<Type>phiprey(wlength); // Available prey
    phiprey = predkernelPP*phiPP+predkernel*phiN;

    array<Type>f(nspecies,wlength); // Feeding level

    // Calculate the feeding level

    for(int k=0;k<nspecies;k++){ // Loop over species
            for(int j=0;j<wlength;j++){ // Calculate the available N
            f(k,j) = Type(1.0)-(IntakeMax(k,j)/SearchVol(k,j)/(phiprey(j)+(IntakeMax(k,j)/SearchVol(k,j))));
            }
    }

    // Calculate mortality
    array<Type> Mtwo(nspecies ,wlength);
    vector<Type>MtwoPP(wPPlength);
    matrix<Type>fN(nspecies,wlength);

    for(int k=0;k<nspecies;k++){ // Loop over species
      for(int j=0;j<wlength;j++){ // Calculate the available N)
        fN(k,j) = N(k,j)*(1.0-f(k,j))*dw(j)*SearchVol(k,j);
      }
    }

    matrix<Type> M2tmp = fN*predkernel;
    matrix<Type> MPP = fN*predkernelPP;
    array<Type> M2(nspecies,wlength);
    array<Type>energy(nspecies,wlength);
    array<Type>SSBm(nspecies,wlength);
    array<Type>gg(nspecies,wlength);
    array<Type>Z(nspecies,wlength);

    // Population matrices for Von Voerster
    matrix<Type>A(nspecies,wlength);
    matrix<Type>B(nspecies,wlength);
    matrix<Type>S(nspecies,wlength);


    for(int k=0;k<nspecies;k++){ // Loop over species


        for(int j=0;j<wlength;j++){ // Loop over lengths
              vector<Type> m1_row1 = M2tmp.col(j);    // 1st row of m1
              M2(k,j) = v(k)*m1_row1.sum(); // Sum all contributions to mortality
              energy(k,j) = alpha(k)*f(k,j)*IntakeMax(k,j)-stdmetab(k,j)-activity(k,j); // Available energy
              SSBm(k,j) = psi(k,j)*energy(k,j); // Spawning rule
              gg(k,j) = energy(k,j)-SSBm(k,j); // Remaining energy used for growth
              Z(k,j) = Z0(k) + M2(k,j)+ Fin(k,j);
        }
// Set up matrices
    for(int idx=1;idx<wlength;idx++){ // Calculate the available N
    A(k,idx) = -Type(1.0)*gg(k,(idx-1))*dt/dw(idx);
    B(k,idx) = 1 + gg(k,idx)*dt/dw(idx) + Z(idx)*dt;
    S(k,idx) = N(k,idx);
    }

    A(k,0) = 0;
    B(k,0) = 0;
    S(k,0) = 0;

REPORT(A)
REPORT(B)
REPORT(S)

    }
//Egg <- sum(param$eRepro[i]*SSBm*N[i,]*dw)   # Egg production (mass/time)
//Rp <- 0.5*param$rho[i]*Egg/w[1]# * param$R0[i]*exp(rnorm(n = 1, mean = 0, sd = param$R.sd[i])) # Egg flux (numbers/time)
//
//# switch param$nRecruitmentType Maybe fix if needed
//# case 0,     # Physiological
//# R <- Rp;    # recruitment <- egg production
//# B(i,1) <- 1 + gg(i,1)$*dt$/dw(1) + Z(1)*dt;
//# N(i,1) <- (N(i,1)+R*dt/dw(1))/B(i,1);
//# case 1,     # Fixed
//# R <- param$fN0(i)*(gg(i,1)*w(1));
//# N(i,1) <- param$fN0(i);
//
//# Beverton-Holt
//R <- param$Rmax[i]*Rp/(param$Rmax[i]+Rp)
//B[i,1] <- 1 + gg[i,1]*dt/dw[1] + Z[1]*dt
//N[i,1] <- (N[i,1]+R*dt/dw[1])/B[i,1]
//
//#
//# Invert matrix
//#
//for (j in 2:nGrid){
//N[i,j] <- (S[i,j]-A[i,j]*N[i,j-1])/B[i,j]
//}



//    Type Btemp = 0.0;
//    Type Ctemp = 0.0;
//    Type SSBtemp = 0.0;
//     for(int j=0;j<wlength;j++){
//        Btemp += N(j)*w(j)*dw(j);
//        Ctemp += Fin(j,i)*N(j)*w(j)*dw(j);
//        SSBtemp += phiMat(j)*N(j)*w(j)*dw(j);
//        Nsave(j,i) = N(j);
//     }
//
//     Catchest(i) = Ctemp;
//     Bio(i) = Btemp;
//     SSB(i) = SSBtemp;


//}  // End time loop
//
//
// Calculate stuff
//vector<Type> Ntemp(nlength);
//vector<Type> Nest(nobs);
//
//for(int k=0;k<(iTimeMax);k++){ // bin the data
//
//   for(int j=0;j<(nlength);j++){
//        Ntemp(j) = 0.0;
////
//        for(int i=0;i<(wlength);i++){
////
//        Type leftEdge = mbins(j);
//        Type rightEdge = mbins(j+1);
////
//       if (w(i) >= leftEdge && w(i)<= rightEdge){
//            Ntemp(j) += Nsave(i,k)*surveysel(i)*catchability; // Fix catchability later
////
//        }
//
//    }
//}
//
//int idxtmp;
//for(int u=0;u<(nlength);u++){
//    idxtmp = yearidx(k)+u;
//    Nest(idxtmp) = Ntemp(u);
//
//}
////
//}

//using namespace density;
//
//  for(int i=0;i<iTimeMax;i++){
//     ans+= -dnorm(logF(i),Type(0.0),SDF,TRUE); // F-Process likelihood
//  }
////
// for(int i=0;i<iTimeMax;i++){
//     ans+= -dnorm(logM(i),Type(0.0), SDM, TRUE); // M-Process likelihood
//  }
//
////   for(int i=1;i<iTimeMax;i++){
////     ans+= -dnorm(logM(i),logM(i-1), SDM, TRUE); // M-Process likelihood
////  }
//
//
// vector<Type> RBH(iTimeMax); //  stock recruitment
//
//  for(int i=0;i<iTimeMax;i++){
//     RBH(i) = (Rmax*Rp(i))/(Rmax+Rp(i));
//     ans+= -dnorm(log(R(i)),log(RBH(i)),SDR,TRUE); // R-Process likelihood
//  }
//// for(int i=1;i<iTimeMax;i++){
////     ans+= -dnorm(logR(i),Type(0.0), SDR, TRUE); // M-Process likelihood
////  }
//
//
//
//for (int i=0;i<nobs;i++){
//    ans += -dnorm(log(Nest(i)+kappa), log(survey(i)+kappa), SDsurv, TRUE); // lognormal fit
//    }
//
//
//for (int i=0;i<iTimeMax;i++){
//    ans += -dnorm(log(Catchest(i)+kappa), log(Catch(i)+kappa), SD, TRUE); // lognormal fit
//    }
//
//// Report calculations
//ADREPORT(Catchest)
//ADREPORT(Bio)
//ADREPORT(Fsave)
//ADREPORT(Msave)
//ADREPORT(R)
//ADREPORT(Nsave)
//ADREPORT(Nest)
//ADREPORT(Fsel)
//ADREPORT(SSB)
//ADREPORT(surveysel)
//ADREPORT(logM)
 // Calculate the likelihood
Type ans = 0.0;
ans = ans + (1-a)*(1-a);

  return 0;
}





//# B(i,1) <- 1 + gg(i,1)$*dt$/dw(1) + Z(1)*dt;
//# N(i,1) <- (N(i,1)+R*dt/dw(1))/B(i,1);
//# case 1,     # Fixed
//# R <- param$fN0(i)*(gg(i,1)*w(1));
//# N(i,1) <- param$fN0(i);
//
//# Beverton-Holt
//R <- param$Rmax[i]*Rp/(param$Rmax[i]+Rp)
//B[i,1] <- 1 + gg[i,1]*dt/dw[1] + Z[1]*dt
//N[i,1] <- (N[i,1]+R*dt/dw[1])/B[i,1]
//
//#
//# Invert matrix
//#
//for (j in 2:nGrid){
//N[i,j] <- (S[i,j]-A[i,j]*N[i,j-1])/B[i,j]
//}



//    Type Btemp = 0.0;
//    Type Ctemp = 0.0;
//    Type SSBtemp = 0.0;
//     for(int j=0;j<wlength;j++){
//        Btemp += N(j)*w(j)*dw(j);
//        Ctemp += Fin(j,i)*N(j)*w(j)*dw(j);
//        SSBtemp += phiMat(j)*N(j)*w(j)*dw(j);
//        Nsave(j,i) = N(j);
//     }
//
//     Catchest(i) = Ctemp;
//     Bio(i) = Btemp;
//     SSB(i) = SSBtemp;


//}  // End time loop
//
//
// Calculate stuff
//vector<Type> Ntemp(nlength);
//vector<Type> Nest(nobs);
//
//for(int k=0;k<(iTimeMax);k++){ // bin the data
//
//   for(int j=0;j<(nlength);j++){
//        Ntemp(j) = 0.0;
////
//        for(int i=0;i<(wlength);i++){
////
//        Type leftEdge = mbins(j);
//        Type rightEdge = mbins(j+1);
////
//       if (w(i) >= leftEdge && w(i)<= rightEdge){
//            Ntemp(j) += Nsave(i,k)*surveysel(i)*catchability; // Fix catchability later
////
//        }
//
//    }
//}
//
//int idxtmp;
//for(int u=0;u<(nlength);u++){
//    idxtmp = yearidx(k)+u;
//    Nest(idxtmp) = Ntemp(u);
//
//}
////
//}

//using namespace density;
//
//  for(int i=0;i<iTimeMax;i++){
//     ans+= -dnorm(logF(i),Type(0.0),SDF,TRUE); // F-Process likelihood
//  }
////
// for(int i=0;i<iTimeMax;i++){
//     ans+= -dnorm(logM(i),Type(0.0), SDM, TRUE); // M-Process likelihood
//  }
//
////   for(int i=1;i<iTimeMax;i++){
////     ans+= -dnorm(logM(i),logM(i-1), SDM, TRUE); // M-Process likelihood
////  }
//
//
// vector<Type> RBH(iTimeMax); //  stock recruitment
//
//  for(int i=0;i<iTimeMax;i++){
//     RBH(i) = (Rmax*Rp(i))/(Rmax+Rp(i));
//     ans+= -dnorm(log(R(i)),log(RBH(i)),SDR,TRUE); // R-Process likelihood
//  }
//// for(int i=1;i<iTimeMax;i++){
////     ans+= -dnorm(logR(i),Type(0.0), SDR, TRUE); // M-Process likelihood
////  }
//
//
//
//for (int i=0;i<nobs;i++){
//    ans += -dnorm(log(Nest(i)+kappa), log(survey(i)+kappa), SDsurv, TRUE); // lognormal fit
//    }
//
//
//for (int i=0;i<iTimeMax;i++){
//    ans += -dnorm(log(Catchest(i)+kappa), log(Catch(i)+kappa), SD, TRUE); // lognormal fit
//    }
//
//// Report calculations
//ADREPORT(Catchest)
//ADREPORT(Bio)
//ADREPORT(Fsave)
//ADREPORT(Msave)
//ADREPORT(R)
//ADREPORT(Nsave)
//ADREPORT(Nest)
//ADREPORT(Fsel)
//ADREPORT(SSB)
//ADREPORT(surveysel)
//ADREPORT(logM)
 // Calculate the likelihood
Type ans = 0.0;
ans = ans + (1-a)*(1-a);

  return 0;
}




