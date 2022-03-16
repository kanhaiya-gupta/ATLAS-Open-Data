#include "GetGaussLimit.C"
#include "GetExpectedLimit.C"
#include "GetBayesBgrCL.C"
#include "GetCLs.C"

Double_t GetBgrCL(Int_t iType,Double_t sig,Double_t bgr,
		  Int_t nobs,Double_t *dcl,TRandom *rnd) {
  Double_t cl;
  if(iType==0) {
    cl = GetBayesBgrCL(sig,bgr,nobs,dcl,rnd);
  } else {
    cl = 1.0-GetCLs(sig,bgr,nobs);
  }
  return cl;
}

void GetLimitBgr(Double_t cl,Double_t bgr,Int_t n0) {
  Double_t sigGauss1=GetGaussLimit(cl,n0+1)-bgr;
  Double_t sig[2],p[2],dp[2];
  TRandom rnd;
  for(Int_t i=0;i<2;i++) {
    Double_t sig0=n0-bgr;
    if(sig0<0.0) sig0=0.;
    Double_t cl0=GetBgrCL(i,sig0,bgr,n0,0,&rnd);
    Double_t sig1=sigGauss1;
    while(sig1<=sig0) sig1+=TMath::Sqrt(n0+0.5);
    Double_t cl1=GetBgrCL(i,sig1,bgr,n0,0,&rnd);
    Double_t clx,dclx;
    do {
      sig[i]=sig0+(sig1-sig0)/(cl1-cl0)*(cl-cl0);
      dclx=0.001*cl;
      clx=GetBgrCL(i,sig[i],bgr,n0,&dclx,&rnd);
      if(TMath::Abs(cl0-clx)<TMath::Abs(cl1-clx)) {
	cl1=clx;
	sig1=sig[i];
      } else {
	cl0=clx;
	sig0=sig[i];
      }
    }
    while(TMath::Abs(clx-cl)>dclx);
    if(i==0) {
      cout<<"Bayes: "<<sig[i]<<"\n";
    } else {
      cout<<"CLs: "<<sig[i]<<"\n";
    }
  }
  cout<<"Poisson: "<<GetPoissonLimit(cl,n0)-bgr<<"\n";
  cout<<"Poisson(expected): "<<GetExpectedLimit(cl,bgr)-bgr<<"\n";
  cout<<"Gaussian approx: "<<GetGaussLimit(cl,n0)-bgr<<"\n";
}
