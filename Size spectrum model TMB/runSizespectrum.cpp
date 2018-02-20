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
  DATA_INTEGER(iTimeMax);

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
  DATA_VECTOR(eRepro); // Recruitment efficiency
  DATA_VECTOR(Rmax); // Carrying capacity
  DATA_VECTOR(Fzero);

//  // Array data
  DATA_ARRAY(Sequ); // Initial size distribution
  DATA_ARRAY(psi); // Mortality ogive
  DATA_MATRIX(predkernel); // Predation kernels per species
  DATA_MATRIX(predkernelPP); // Predation kernel for PP
//
//  // Scalar data
//  // DATA_SCALAR(Zbase);
  DATA_SCALAR(n);
  DATA_SCALAR(p);
  DATA_SCALAR(q); // Search volume exponent
  DATA_SCALAR(dt);
  DATA_SCALAR(rR); // Growth rate of PP
  DATA_SCALAR(lR); // Metabolic constant of PP
  DATA_SCALAR(kR); // slope of PP
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


array<Type> Fin(nspecies,wlength); // Survey selectivity
////
for(int k=0;k<nspecies;k++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Initialize vectors
            // Fishing selectivity
        Fin(k,i) =Fzero(k) * pow( (Type(1.0)+pow( (w(i)/(nF(k) * wInf(k))),-muF(k))),Type(-1));
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


vector<Type> nPP0(wPPlength); // Background spectrum
vector<Type> rrPP(wPPlength); // Growth rate at size
 for(int i=0;i<wPPlength;i++){ // Loop over PP sizes

  nPP0(i) = kappaR * pow(wPP(i),kR);
  rrPP(i) = rR * pow(wPP(i),-lR); // Weight specific growth rate

}

vector<Type>nPP = nPP0;


array<Type> N(nspecies,wlength); // Size spectrum
for(int k=0;k<nspecies;k++){ // Loop over species
        for(int i=0;i<wlength;i++){ // Initialize vectors
                N(k,i) = Sequ(k,i);
                }
            }

// things to save
array<Type>R(nspecies,iTimeMax);
array<Type>Rpsave(nspecies,iTimeMax);

array<Type>Nsave(nspecies,wlength,iTimeMax);
array<Type>Csave(nspecies,iTimeMax);
array<Type>Bsave(npecies,iTimemax);

//int i = 0;
for(int i=0;i<iTimeMax;i++){ // start time loop

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
    matrix<Type> MPPtmp = fN*predkernelPP;
    vector<Type>MPP(wPPlength);
    for(int j=0;j<wPPlength;j++){ // Calculate the available PP
              vector<Type> m2_col2 = MPPtmp.col(j);    // 1st row of m1
              MPP(j) = m2_col2.sum(); // Sum all contributions to mortality
        }

    array<Type> M2(nspecies,wlength);

    array<Type>energy(nspecies,wlength);
    array<Type>SSBm(nspecies,wlength);
    array<Type>gg(nspecies,wlength);
    array<Type>Ztot(nspecies,wlength);


    // Population matrices for Von Voerster
    array<Type>Amat(nspecies,wlength);
    array<Type>B(nspecies,wlength);
    array<Type>S(nspecies,wlength);

    for(int k=0;k<nspecies;k++){ // Loop over species


        for(int j=0;j<wlength;j++){ // Loop over lengths
              vector<Type> m1_row1 = M2tmp.col(j);    // 1st row of m1
              M2(k,j) = v(k)*m1_row1.sum(); // Sum all contributions to mortality
              energy(k,j) = alpha(k)*f(k,j)*IntakeMax(k,j)-stdmetab(k,j)-activity(k,j); // Available energy
              if(energy(k,j)<0){
                energy(k,j) = 0.0;
              }

              SSBm(k,j) = psi(k,j)*energy(k,j); // Spawning rule
              gg(k,j) = energy(k,j)-SSBm(k,j); // Remaining energy used for growth
              Ztot(k,j) = Z0(k) + M2(k,j)+ Fin(k,j);//
        }

// Set up matrices
    for(int j=1;j<wlength;j++){ // Calculate the available N
     Amat(k,j) = -gg(k,j-1)* dt/dw(j);
     B(k,j) = Type(1.0) + gg(k,j)*dt/dw(j) + Ztot(k,j)*dt;
     S(k,j) = N(k,j);
    }

    Type Egg = 0.0;
    for(int j=0;j<wlength;j++){ // Loop over lengths
    Egg += eRepro(k)*SSBm(k,j)*N(k,j)*dw(j);   // Egg production (mass/time)
    }
    Type Rp = 0.5*Egg/w(0);

    Type Rtemp = (Rmax(k)*Rp/(Rmax(k)+Rp)); //;
    R(k,i) = Rtemp;
    Rpsave(k,i)= Rp;


    B(k,0) = Type(1.0) + gg(k,0)*dt/dw(0) + Ztot(k,0)*dt;
    N(k,0) = (N(k,0)+ Rtemp*dt/dw(0))/B(0);

    for(int j=1;j<wlength;j++){ // Invert matrices
      N(k,j) = (S(k,j)-Amat(k,j)*N(k,j-1))/B(k,j);

    }
    Type Btemp = 0.0;
    Type Ctemp = 0.0;
    Type SSBtemp = 0.0;
     for(int j=0;j<wlength;j++){
        Bsave(k,i) += N(j)*w(j)*dw(j);
        Csave(k,i) += Fin(k,j)*N(k,j)*w(j)*dw(j);
        Nsave(k,j,i) = N(k,j);
     }


    }

    vector<Type> tmp(wPPlength); // Correct primary production
    for(int j=0;j<wPPlength;j++){ // Loop over PP sizes
        tmp(j) = (rrPP(j)*nPP0(j) / (rrPP(j) + MPP(j)));
        nPP(j) = tmp(j) - (tmp(j) - nPP(j))*exp(-(rrPP(j)+MPP(j))*dt);
    }

    // Save some things




} // End of time loop




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

Type ans = 0.0;
ans = ans + (1-a)*(1-a);

REPORT(Nsave)
REPORT(R)
REPORT(Rpsave)

  return 0;
}




