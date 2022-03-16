#include "GetBayesCL.C"
#include "GetPoissonLimit.C"
#include "GetGaussLimit.C"

void GetLimit(Double_t cl,Int_t n0) {
  Double_t muGauss1=GetGaussLimit(cl,n0+1);
  Double_t muBayes;
  TRandom3 rnd;
  Double_t mu0=n0;
  Double_t cl0=GetBayesCL(mu0,n0,0,&rnd);
  Double_t mu1=muGauss1;
  Double_t cl1=GetBayesCL(mu1,n0,0,&rnd);
  Double_t clx,dclx;
  do {
    muBayes=mu0+(mu1-mu0)/(cl1-cl0)*(cl-cl0);
    clx=GetBayesCL(muBayes,n0,&dclx,&rnd);
    if(TMath::Abs(cl0-clx)<TMath::Abs(cl1-clx)) {
      cl1=clx;
      mu1=muBayes;
    } else {
      cl0=clx;
      mu0=muBayes;
    }
  }
  while(TMath::Abs(clx-cl)>dclx);
  cout<<"Bayesian limit: "<<muBayes<<"\n";
  cout<<"Frequentist limit: "<<GetPoissonLimit(cl,n0)<<"\n";
  cout<<"Gaussian approx: "<<GetGaussLimit(cl,n0)<<"\n";
}
