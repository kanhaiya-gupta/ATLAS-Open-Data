Double_t GetPoissonCL(Double_t mu,Int_t n0,Double_t *dp=0,TRandom *rnd=0) {
  /* calculates the Poisson probability (1-CL) to observe N>n0
     given mu. It uses a Monte Carlo technique.
     If dp!=0, also return the error */
  Int_t nTry=100000;
  Double_t nSuccess=0;
  TRandom *ownedRnd=0;
  if(!rnd) ownedRnd=new TRandom3();
  if(!rnd) rnd=ownedRnd;
  for(Int_t iTry=0;iTry<nTry;iTry++) {
    Int_t nObs=rnd->Poisson(mu);
    if(nObs>n0) {
      nSuccess += 1.0;
    }
  }
  Double_t p=nSuccess/(Double_t)nTry;
  if(dp) *dp=TMath::Sqrt((nSuccess)*(nTry-nSuccess))/TMath::Power(nTry,1.5);
  if(ownedRnd) delete ownedRnd;
  return p;
}
