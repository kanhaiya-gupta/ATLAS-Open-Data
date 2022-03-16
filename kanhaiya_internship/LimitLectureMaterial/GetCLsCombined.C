Double_t GetX(Int_t nChan,Int_t *nobs,Double_t signal,Double_t *eff,Double_t *b) {
  Double_t x=0.;
  for(Int_t i=0;i<nChan;i++) {
    Double_t s=signal*eff[i];
    x += nobs[i]*s/b[i];
  }
  return x;
}

void GetCLsCombined(Double_t signal,TRandom *rnd=0) {
  TRandom *ownedRnd=0;
  if(!rnd) ownedRnd=new TRandom3();
  if(!rnd) rnd=ownedRnd;
  // data to be tested
  static Int_t nObs[2]  ={7   , 2  };
  static Double_t bgr[2]={6.5 ,1.8 };
  static Double_t eff[2]={0.5 , 0.5};
  // number of toy experiments
  Int_t nTry=100000;
  // count toy wrt data experiments
  Double_t ndata_sb=0.0;
  Double_t ndata_b=0.0;
  // count toy wrt toy experiments for "expected" CL
  Double_t nexpect_sb=0.0;
  Double_t nexpect_b=0.0;
  // observed X from data
  Double_t Xobs=GetX(2,nObs,signal,eff,bgr);
  // toy experiments
  for(Int_t iTry=0;iTry<nTry;iTry++) {
    Int_t n_b[2],n_sb[2],n_b2[2];
    for(Int_t i=0;i<2;i++) {
      n_b[i]=rnd->Poisson(bgr[i]);
      n_b2[i]=rnd->Poisson(bgr[i]);
      n_sb[i]=rnd->Poisson(signal*eff[i]+bgr[i]);
    }
    Double_t X_b=GetX(2,n_b,signal  ,eff,bgr);
    Double_t X_sb=GetX(2,n_sb,signal,eff,bgr);
    Double_t Xtoy=GetX(2,n_b2,signal,eff,bgr);
    if(X_b<=Xobs) ndata_b++;
    if(X_sb<=Xobs) ndata_sb++;
    if(X_b<=Xtoy) nexpect_b++;
    if(X_sb<=Xtoy) nexpect_sb++;
  }
  if(ownedRnd) delete ownedRnd;
  Double_t dataCL_s=ndata_sb/ndata_b;
  Double_t expectCL_s=nexpect_sb/nexpect_b;
  cout<<"CLS(data)="<<dataCL_s<<" CLS(expected)="<<expectCL_s<<"\n";
}
